"""Behavior lifecycle tests; never open hardware or send network traffic."""
import itertools
from pathlib import Path
import queue
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

import behavior_gui as gui


class FakeSerial:
    def __init__(self):
        self.input = bytearray()
        self.written = []
        self.closed = False

    @property
    def in_waiting(self):
        return len(self.input)

    def read(self, count):
        result = bytes(self.input[:count])
        del self.input[:count]
        return result

    def write(self, data):
        self.written.append(data)
        return len(data)

    def close(self):
        self.closed = True


class BehaviorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.service = gui.BehaviorService()
        self.addCleanup(self.service.do_disconnect)

    def module(self, source):
        path = self.base / "event protocol.py"
        path.write_text(source, encoding="utf-8")
        return path

    def sender(self):
        sock = Mock()
        with patch.object(gui.socket, "socket", return_value=sock):
            sender = gui.UdpSender("127.0.0.1", 5005, Mock(), itertools.count(100).__next__)
        self.addCleanup(sender.close)
        return sender, sock

    def logs(self):
        messages = []
        while True:
            try:
                messages.append(self.service.logs.get_nowait())
            except queue.Empty:
                return "\n".join(messages)

    def test_udp_repeats_share_sequence_and_close_cancels_delayed_work(self):
        with patch.object(gui.time, "monotonic", return_value=10) as clock:
            sender, sock = self.sender()
            self.assertEqual(sender.send("CMD PID_ON"), 100)
            sender.call_later("CMD PID_OFF", 0.5)
            clock.return_value = 10.01
            sender.poll()
            self.assertEqual(sock.sendto.call_count, 3)
            self.assertEqual({c.args[0] for c in sock.sendto.call_args_list}, {b"MSG 100 CMD PID_ON\n"})
            sender.close()
            clock.return_value = 11
            sender.poll()
            sender.send("CMD SET_TARGET 1")
            self.assertEqual(sock.sendto.call_count, 3)
            self.assertFalse(sender.pending)

    def test_udp_delayed_command_gets_new_sequence(self):
        with patch.object(gui.time, "monotonic", return_value=10) as clock:
            sender, sock = self.sender()
            sender.call_later("CMD PID_OFF", 0.5)
            sender.poll()
            sock.sendto.assert_not_called()
            clock.return_value = 10.5
            sender.poll()
            sock.sendto.assert_called_once_with(b"MSG 100 CMD PID_OFF\n", ("127.0.0.1", 5005))

    def test_protocol_import_skips_main_and_redirects_legacy_timer(self):
        path = self.module('''from dataclasses import dataclass
@dataclass
class Event:
    type: str
def parse_arduino_line(line):
    return Event(line)
def delay_send(*args):
    raise AssertionError("Legacy timer must be replaced")
def apply_session_logic(event, udp):
    delay_send(udp, event.type, 0.4)
if __name__ == "__main__":
    raise AssertionError("Must not launch a second GUI")
''')
        protocol = gui.Protocol(path, Mock())
        sender = Mock()
        protocol.handle("CMD PID_ON", sender)
        sender.call_later.assert_called_once_with("CMD PID_ON", 0.4)
        protocol.close()
        self.assertNotIn(protocol.name, sys.modules)

    def test_invalid_protocol_leaves_no_imported_module_or_udp_socket(self):
        path = self.module("x = 1\n")
        self.service.ser = FakeSerial()
        before = set(sys.modules)
        with patch.object(gui.socket, "socket") as socket_factory:
            with self.assertRaisesRegex(ValueError, "parse_arduino_line"):
                self.service.do_start_protocol(gui.Settings(protocol=str(path)))
            socket_factory.return_value.close.assert_called_once()
        self.assertFalse([n for n in set(sys.modules) - before if n.startswith("_behavior_protocol_")])
        self.assertIsNone(self.service.protocol)
        self.assertIsNone(self.service.udp)

    def test_protocol_restarts_with_fresh_globals_and_source(self):
        path = self.module('''count = 0
def parse_arduino_line(line): return line
def apply_session_logic(event, udp):
    global count
    count += 1
''')
        first = gui.Protocol(path, Mock())
        first.handle("cue", Mock())
        self.assertEqual(first.module.count, 1)
        first.close()
        second = gui.Protocol(path, Mock())
        self.assertEqual(second.module.count, 0)
        second.close()

    def test_serial_partial_lines_preserve_events_and_raw_text(self):
        self.service.ser = FakeSerial()
        protocol = Mock()
        protocol.handle.return_value = None
        self.service.protocol = protocol
        self.service.ser.input.extend(b"Trial: 1 Cue st")
        self.service._read_serial()
        protocol.handle.assert_not_called()
        self.service.ser.input.extend(b"art (Pair #1) Time: 1.0\r\nunknown\n")
        self.service._read_serial()
        self.assertEqual([c.args[0] for c in protocol.handle.call_args_list],
                         ["Trial: 1 Cue start (Pair #1) Time: 1.0", "unknown"])
        self.assertIn("[RX] unknown", self.logs())

    def test_serial_commands_preserve_spaces_and_line_endings(self):
        self.service.ser = FakeSerial()
        self.service.do_send("  command  ", "Both (CRLF)")
        self.service.do_send("x", "None")
        self.assertEqual(self.service.ser.written, [b"  command  \r\n", b"x"])

    def test_failed_serial_write_disconnects_without_retry(self):
        serial = FakeSerial()
        serial.write = Mock(side_effect=OSError("unplugged"))
        self.service.ser = serial
        with self.assertRaises(OSError):
            self.service.do_send("start", "Newline (LF)")
        self.assertTrue(serial.closed)
        self.assertIsNone(self.service.ser)
        serial.write.assert_called_once()

    def test_protocol_error_stops_udp_but_keeps_raw_monitor(self):
        serial = FakeSerial()
        self.service.ser = serial
        self.service.protocol = Mock()
        self.service.protocol.handle.side_effect = ValueError("bad event")
        udp, _ = self.sender()
        udp.call_later("CMD PID_OFF", 10)
        self.service.udp = udp
        serial.input.extend(b"bad line\nnext line\n")
        self.service._read_serial()
        self.assertIsNone(self.service.protocol)
        self.assertIs(self.service.ser, serial)
        self.assertTrue(udp.closed)
        self.assertFalse(udp.pending)
        self.assertIn("[RX] next line", self.logs())

    def sketch_settings(self):
        folder = self.base / "Sketch with spaces"
        folder.mkdir()
        path = folder / (folder.name + ".ino")
        path.write_text("void setup() {}\nvoid loop() {}\n")
        return gui.Settings(sketch=str(path), board="arduino:avr:uno", port="COM4", reconnect=False)

    def test_upload_closes_serial_and_protocol_before_compile(self):
        settings = self.sketch_settings()
        serial = FakeSerial()
        self.service.ser = serial
        self.service.protocol = Mock()
        sender, _ = self.sender()
        sender.call_later("CMD PID_OFF", 1)
        self.service.udp = sender
        commands = []
        def run(args, timeout):
            self.assertTrue(serial.closed)
            self.assertTrue(sender.closed)
            self.assertIsNone(self.service.protocol)
            commands.append(args)
        with patch.object(gui, "cli_executable", return_value="/path with spaces/arduino-cli"), \
                patch.object(self.service, "_run_cli", side_effect=run):
            self.service.do_upload(settings)
        self.assertEqual([c[1] for c in commands], ["compile", "upload"])
        self.assertEqual(commands[0][-1], str(Path(settings.sketch).resolve().parent))
        self.assertEqual(commands[0][commands[0].index("--build-path") + 1],
                         commands[1][commands[1].index("--build-path") + 1])
        self.assertFalse(self.service.busy)

    def test_compile_failure_never_uploads_or_reconnects(self):
        settings = self.sketch_settings()
        settings.reconnect = True
        with patch.object(gui, "cli_executable", return_value="cli"), \
                patch.object(self.service, "_run_cli", side_effect=RuntimeError("compile failed")) as run, \
                patch.object(self.service, "do_connect") as connect:
            with self.assertRaisesRegex(RuntimeError, "compile failed"):
                self.service.do_upload(settings)
        run.assert_called_once()
        connect.assert_not_called()
        self.assertFalse(self.service.busy)

    def test_bad_sketch_does_not_disconnect_existing_session(self):
        self.service.ser = FakeSerial()
        with patch.object(gui, "cli_executable", return_value="cli"):
            with self.assertRaises(ValueError):
                self.service.do_upload(gui.Settings(sketch=str(self.base / "missing.ino")))
        self.assertFalse(self.service.ser.closed)

    def test_cli_process_streams_errors_and_rejects_failure(self):
        script = self.base / "fake cli.py"
        script.write_text("import sys\nprint('missing board core', flush=True)\nsys.exit(2)\n")
        with self.assertRaisesRegex(RuntimeError, "exit 2"):
            self.service._run_cli([sys.executable, str(script)], 5)
        self.assertIn("missing board core", self.logs())

    def test_bundled_rpe_emits_original_cue_reward_and_omission_commands(self):
        protocol = gui.Protocol(gui.Settings().protocol, Mock())
        self.addCleanup(protocol.close)
        udp = Mock()
        event = protocol.handle("Trial: 12 Cue start (Pair #3) Time: 123.45", udp)
        self.assertEqual((event.type, event.trial), ("TRIAL_ON", 12))
        udp.send.assert_called_with("CMD PID_ON 200")
        udp.call_later.assert_called_with("CMD SET_TARGET 0.8000 0.25 natural", 0.4)
        protocol.handle("Trial: 12 Reward block: big reward Time: 124.0", udp)
        udp.send.assert_called_with("CMD SET_TARGET 0.8000")
        udp.call_later.assert_called_with("CMD PID_OFF 500", 0.5)
        protocol.handle("Trial: 13 Reward block: omission Time: 130.0", udp)
        udp.send.assert_called_with("CMD SET_TARGET -0.1000 0.25 natural")
        self.assertIsNone(protocol.handle("unrecognized message", udp))

    def test_service_shutdown_closes_connections(self):
        serial = FakeSerial()
        self.service.ser = serial
        udp, _ = self.sender()
        self.service.udp = udp
        self.service.thread.start()
        self.service.shutdown.set()
        self.service.thread.join(timeout=2)
        self.assertFalse(self.service.thread.is_alive())
        self.assertTrue(serial.closed)
        self.assertTrue(udp.closed)


if __name__ == "__main__":
    unittest.main()
