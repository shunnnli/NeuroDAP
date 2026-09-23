#!/usr/bin/env python3
"""Arduino upload, serial monitor, and optional Python-to-UDP behavior bridge.

Run: python gui_behavior.py
Dependencies/setup and the protocol API: gui/README.md.
Importing this module does not open a window or connect to hardware.
"""
from __future__ import annotations

import argparse
import heapq
import importlib.util
import ipaddress
import json
from pathlib import Path
import queue
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time
import uuid
from dataclasses import asdict, dataclass, field, fields

from .parameters import ParameterEditor, apply_overrides, prepare_sketch

GUI_DIR = Path(__file__).resolve().parent
ROOT = GUI_DIR.parent
PROFILES_DIR = GUI_DIR / "profiles"
PARSER_DIR = ROOT.parent / "BrainClamp" / "scripts"
LINE_ENDINGS = {"Newline (LF)": "\n", "Carriage return (CR)": "\r",
                "Both (CRLF)": "\r\n", "None": ""}
BOARDS = ("arduino:avr:uno", "arduino:avr:mega:cpu=atmega2560",
          "arduino:avr:nano:cpu=atmega328", "arduino:avr:nano:cpu=atmega328old",
          "arduino:avr:leonardo", "arduino:samd:mkrzero")


@dataclass
class Settings:
    sketch: str = ""
    board: str = "arduino:avr:mega:cpu=atmega2560"
    port: str = "COM4"
    baud: str = "115200"
    cli: str = "arduino-cli"
    protocol: str = str(PARSER_DIR / "send_event_RPE.py")
    host: str = "192.168.50.132"
    udp_port: str = "5005"
    line_ending: str = "Newline (LF)"
    reconnect: bool = True
    parameter_overrides: dict[str, str] = field(default_factory=dict)
    parameter_source_hash: str = ""


def serial_module():
    try:
        import serial
    except ImportError as exc:
        raise RuntimeError("pySerial is missing. Install with: python -m pip install pyserial") from exc
    return serial


def cli_executable(value):
    value = str(Path(value).expanduser()) if value else "arduino-cli"
    found = shutil.which(value)
    if not found:
        raise ValueError("Arduino CLI was not found. Install it or select its executable in the GUI.")
    return found


def sketch_directory(value):
    if not value:
        raise ValueError("Select an Arduino .ino file first.")
    path = Path(value).expanduser().resolve()
    if path.is_file() and path.suffix.lower() == ".ino":
        path = path.parent
    if not path.is_dir() or not (path / (path.name + ".ino")).is_file():
        raise ValueError("The sketch folder must contain a main .ino file with the same name as the folder.")
    return path


def serial_settings(settings):
    if not settings.port.strip():
        raise ValueError("Select a serial port.")
    baud = int(settings.baud)
    if baud <= 0:
        raise ValueError("Baud rate must be positive.")
    return settings.port.strip(), baud


class UdpSender:
    """Wire-compatible with send_event_RPE.py; one sequence per logical command.

    Repeats and delayed commands run on the service thread. Closing the sender
    invalidates all scheduled work, including references held by a protocol.
    """
    def __init__(self, host, port, log, next_sequence, resend_n=3, resend_gap_s=0.003):
        ipaddress.IPv4Address(host)
        port = int(port)
        if not 1 <= port <= 65535:
            raise ValueError("UDP port must be between 1 and 65535.")
        self.addr = (host, port)
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        self.sock.setblocking(False)
        self.log = log
        self.next_sequence = next_sequence
        self.resend_n = resend_n
        self.resend_gap_s = resend_gap_s
        self.seq = 0
        self.pending = []
        self.counter = 0
        self.closed = False
        self.lock = threading.RLock()

    def _schedule(self, callback, delay):
        self.counter += 1
        heapq.heappush(self.pending, (time.monotonic() + max(0, float(delay)), self.counter, callback))

    def send(self, payload):
        with self.lock:
            if self.closed:
                return None
            self.seq = self.next_sequence()
            packet = f"MSG {self.seq} {payload}\n".encode("ascii")
            self.sock.sendto(packet, self.addr)
            for n in range(1, self.resend_n):
                self._schedule(lambda p=packet: self.sock.sendto(p, self.addr), n * self.resend_gap_s)
            self.log(f"[UDP sent #{self.seq}] {payload}")
            return self.seq

    def call_later(self, payload, delay_s):
        with self.lock:
            if not self.closed:
                self._schedule(lambda: self.send(payload), delay_s)

    def poll(self):
        with self.lock:
            while not self.closed and self.pending and self.pending[0][0] <= time.monotonic():
                _, _, callback = heapq.heappop(self.pending)
                callback()

    def close(self):
        with self.lock:
            self.closed = True
            self.pending.clear()
            self.sock.close()


class Protocol:
    """Load a fresh module per session, without running its __main__ entry point."""
    def __init__(self, path, log):
        path = Path(path).expanduser().resolve()
        if not path.is_file() or path.suffix != ".py":
            raise ValueError("Select an event-parser Python (.py) file.")
        self.name = "_behavior_protocol_" + uuid.uuid4().hex
        spec = importlib.util.spec_from_file_location(self.name, path)
        self.module = importlib.util.module_from_spec(spec)
        sys.modules[self.name] = self.module  # dataclasses resolve their defining module
        try:
            # Compile source directly so a rapid edit/restart cannot reuse stale bytecode.
            exec(compile(path.read_bytes(), str(path), "exec"), self.module.__dict__)
            for name in ("parse_arduino_line", "apply_session_logic"):
                if not callable(getattr(self.module, name, None)):
                    raise ValueError(f"Protocol must define {name}(). See gui/README.md.")
            # Compatibility with the supplied send_event_RPE.py. Its GUI/main is
            # never called, and its Timer helper is replaced with cancellable work.
            self.module.delay_send = lambda udp, cmd, delay_s: udp.call_later(cmd, delay_s) if udp else None
            if callable(getattr(self.module, "register_serial_monitor_log", None)):
                self.module.register_serial_monitor_log(log)
        except Exception:
            sys.modules.pop(self.name, None)
            raise

    def handle(self, line, udp):
        event = self.module.parse_arduino_line(line)
        if event is not None:
            self.module.apply_session_logic(event, udp)
        return event

    def close(self):
        try:
            if callable(getattr(self.module, "register_serial_monitor_log", None)):
                self.module.register_serial_monitor_log(None)
        finally:
            sys.modules.pop(self.name, None)


class BehaviorService:
    """Own serial/UDP on one worker; communicate with Tk only through queues."""
    def __init__(self):
        self.commands = queue.Queue()
        self.events = queue.Queue()
        self.logs = queue.Queue(maxsize=10000)
        self.shutdown = threading.Event()
        self.ser = None
        self.protocol = None
        self.udp = None
        self.buffer = bytearray()
        self.busy = False
        self.sequence = time.time_ns() // 1_000_000  # survives protocol restarts; receiver accepts integer IDs
        self.dropped = 0
        self.thread = threading.Thread(target=self._run, daemon=True, name="behavior-service")

    def log(self, message):
        try:
            self.logs.put_nowait(f"{time.strftime('%H:%M:%S')} {message}")
        except queue.Full:
            self.dropped += 1

    def state(self):
        self.events.put(("state", {"connected": self.ser is not None,
                                  "protocol": self.protocol is not None, "busy": self.busy}))

    def submit(self, name, *args):
        self.commands.put((name, args))

    def _next_sequence(self):
        self.sequence += 1
        return self.sequence

    def _run(self):
        self.state()
        try:
            while not self.shutdown.is_set():
                try:
                    name, args = self.commands.get(timeout=0.002 if self.ser else 0.05)
                    success = False
                    try:
                        getattr(self, "do_" + name)(*args)
                        success = True
                    except Exception as exc:
                        self.log(f"[error] {exc}")
                        self.events.put(("error", str(exc)))
                    finally:
                        self.state()
                        self.events.put(("done", (name, success)))
                except queue.Empty:
                    pass
                if self.shutdown.is_set():
                    break
                if self.ser:
                    try:
                        self._read_serial()
                    except Exception as exc:
                        self.log(f"[serial error] {exc}")
                        self.do_disconnect()
                        self.state()
                if self.udp:
                    try:
                        self.udp.poll()
                    except Exception as exc:
                        self.log(f"[protocol stopped] {exc}")
                        self.do_stop_protocol()
                        self.state()
        finally:
            self.do_disconnect()

    def do_refresh(self):
        serial_module()
        from serial.tools import list_ports
        ports = list(list_ports.comports())
        self.events.put(("ports", [p.device for p in ports]))
        self.log("[ports] " + (", ".join(f"{p.device}: {p.description}" for p in ports) or "No serial ports found."))

    def do_detect(self, settings):
        cli = cli_executable(settings.cli)
        result = subprocess.run([cli, "board", "list", "--json"], capture_output=True,
                                text=True, encoding="utf-8", errors="replace", timeout=20)
        if result.returncode:
            raise RuntimeError(result.stderr.strip() or result.stdout.strip() or "Board detection failed.")
        data = json.loads(result.stdout)
        ports = data.get("detected_ports", []) if isinstance(data, dict) else data
        matches = []
        for item in ports:
            if item.get("port", {}).get("address") == settings.port:
                matches.extend(b["fqbn"] for b in item.get("matching_boards", []) if b.get("fqbn"))
        if len(matches) == 1:
            self.events.put(("board", matches[0]))
            self.log(f"[board] {matches[0]}")
        else:
            self.log("[board] Could not identify one board. Choose or enter its board ID manually.")

    def do_connect(self, settings):
        port, baud = serial_settings(settings)
        if self.ser:
            raise RuntimeError("Disconnect before opening another serial port.")
        self.ser = serial_module().Serial(port, baud, timeout=0, write_timeout=1)
        self.buffer.clear()
        self.log(f"[connected] {port} @ {baud}. Opening the port may reset the Arduino.")

    def do_stop_protocol(self):
        udp, protocol = self.udp, self.protocol
        self.udp = self.protocol = None
        if udp:
            udp.close()
        if protocol:
            try:
                protocol.close()
            except Exception as exc:
                self.log(f"[protocol cleanup] {exc}")
            self.log("[protocol stopped] Pending UDP commands cancelled. Remote device state is unchanged.")

    def do_disconnect(self):
        self.do_stop_protocol()
        ser, self.ser = self.ser, None
        if ser:
            try:
                ser.close()
            except Exception as exc:
                self.log(f"[serial close] {exc}")
            self.log("[disconnected]")
        self.buffer.clear()

    def do_start_protocol(self, settings):
        if not self.ser:
            raise RuntimeError("Connect the Arduino before starting the protocol.")
        if self.protocol:
            raise RuntimeError("Stop the current protocol first.")
        udp = UdpSender(settings.host.strip(), settings.udp_port, self.log, self._next_sequence)
        try:
            protocol = Protocol(settings.protocol, self.log)
        except Exception:
            udp.close()
            raise
        self.udp, self.protocol = udp, protocol
        self.log(f"[protocol started] {Path(settings.protocol).name} -> {udp.addr[0]}:{udp.addr[1]}")

    def do_send(self, text, ending):
        if not self.ser:
            raise RuntimeError("Connect the Arduino before sending a command.")
        packet = (text + LINE_ENDINGS[ending]).encode("utf-8")
        try:
            count = self.ser.write(packet)
            if count != len(packet):
                raise IOError("Incomplete serial write; command was not retried.")
        except Exception:
            self.do_disconnect()
            raise
        self.log(f"[TX] {text!r} ({ending})")

    def _read_serial(self):
        self.buffer.extend(self.ser.read(min(self.ser.in_waiting, 65536)))
        # Keep incomplete lines across reads; readline(timeout=...) can split events.
        for _ in range(500):
            boundary = self.buffer.find(b"\n")
            if boundary < 0:
                break
            raw = bytes(self.buffer[:boundary]).rstrip(b"\r")
            del self.buffer[:boundary + 1]
            line = raw.decode("utf-8", errors="replace")
            self.log(f"[RX] {line}")
            if self.protocol:
                try:
                    event = self.protocol.handle(line, self.udp)
                    if event is not None:
                        kind = getattr(event, "type", str(event))
                        self.log(f"[event] {kind}")
                except Exception as exc:
                    self.log(f"[protocol stopped] {type(exc).__name__}: {exc}")
                    self.do_stop_protocol()
                    self.state()
        if len(self.buffer) > 1_048_576:
            self.buffer.clear()
            raise RuntimeError("Serial input exceeded 1 MB without enough newline delimiters.")

    def _run_cli(self, args, timeout):
        self.log("[Arduino CLI] " + subprocess.list2cmdline(args))
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, encoding="utf-8", errors="replace")
        def output():
            for line in proc.stdout:
                self.log("[upload] " + line.rstrip())
        reader = threading.Thread(target=output, daemon=True)
        reader.start()
        deadline = time.monotonic() + timeout
        try:
            while proc.poll() is None:
                if self.shutdown.wait(0.05) or time.monotonic() > deadline:
                    proc.terminate()
                    try:
                        proc.wait(timeout=3)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                        proc.wait()
                    raise RuntimeError("Arduino CLI stopped or timed out. Check the board before reconnecting.")
            reader.join(timeout=2)
            if proc.returncode:
                raise RuntimeError(f"Arduino CLI {args[1]} failed (exit {proc.returncode}). See monitor output.")
        finally:
            if not reader.is_alive():
                proc.stdout.close()

    def do_upload(self, settings):
        # Validate before releasing a working serial connection.
        cli = cli_executable(settings.cli)
        sketch = sketch_directory(settings.sketch)
        serial_settings(settings)
        if len(settings.board.strip().split(":")) < 3:
            raise ValueError("Choose a board ID, for example arduino:avr:uno.")
        try:
            with tempfile.TemporaryDirectory(prefix="behavior-build-") as workspace:
                staged = prepare_sketch(sketch, Path(workspace) / "source",
                                        settings.parameter_overrides, settings.parameter_source_hash)
                # Copy and validate first, so invalid/stale edits leave a working
                # serial connection and protocol untouched.
                self.busy = True
                self.do_disconnect()
                self.state()
                self.log(f"[sketch copy] {sketch} -> {staged}")
                for name, value in settings.parameter_overrides.items():
                    self.log(f"[parameter] {name} = {value}")
                build = str(Path(workspace) / "build")
                common = ["--fqbn", settings.board.strip(), "--build-path", build]
                self._run_cli([cli, "compile", *common, str(staged)], 600)
                self._run_cli([cli, "upload", *common, "--port", settings.port.strip(), str(staged)], 120)
            self.log("[upload complete] Protocol is stopped; start it when the board is ready.")
            if settings.reconnect and not self.shutdown.is_set():
                self.shutdown.wait(1)
                deadline = time.monotonic() + 8
                while not self.shutdown.is_set():
                    try:
                        self.do_connect(settings)
                        break
                    except Exception as exc:
                        if time.monotonic() >= deadline:
                            self.log(f"[reconnect failed] {exc}. Refresh ports and connect manually.")
                            break
                        self.shutdown.wait(0.4)
        finally:
            self.busy = False


class BehaviorGUI:
    def __init__(self, root):
        import tkinter as tk
        from tkinter import ttk
        from tkinter.scrolledtext import ScrolledText
        self.tk, self.ttk, self.root = tk, ttk, root
        root.title("NeuroDAP | Behavior control")
        root.geometry("1480x860")
        root.minsize(1230, 780)
        self.service = BehaviorService()
        self.connected = self.protocol_running = self.busy = self.pending = False
        self.closing = False
        self.parameter_sketch = None
        self.parameter_reload_job = None
        self.vars = {f.name: (tk.BooleanVar(value=getattr(Settings(), f.name)) if f.type == "bool"
                             else tk.StringVar(value=getattr(Settings(), f.name))) for f in fields(Settings)
                     if f.name not in ("parameter_overrides", "parameter_source_hash")}
        self.status = tk.StringVar(value="Disconnected")
        self.command = tk.StringVar()
        self.autoscroll = tk.BooleanVar(value=True)
        self.config_widgets, self.protocol_widgets = [], []
        style = ttk.Style(root)
        style.configure("Title.TLabel", font=("Helvetica", 20, "bold"))
        style.configure("Hint.TLabel", foreground="#536477")

        outer = ttk.Frame(root, padding=18)
        outer.grid(sticky="nsew")
        root.rowconfigure(0, weight=1)
        root.columnconfigure(0, weight=1)
        outer.columnconfigure(0, weight=1)
        outer.rowconfigure(4, weight=1)
        header = ttk.Frame(outer)
        header.grid(row=0, column=0, sticky="ew", pady=(0, 16))
        header.columnconfigure(0, weight=1)
        ttk.Label(header, text="Behavior control", style="Title.TLabel").grid(row=0, column=0, sticky="w")
        self.load_button = ttk.Button(header, text="Load profile", command=self.load_profile)
        self.load_button.grid(row=0, column=1, padx=4)
        ttk.Button(header, text="Save profile", command=self.save_profile).grid(row=0, column=2)

        columns = ttk.Frame(outer)
        columns.grid(row=1, column=0, sticky="nsew")
        columns.columnconfigure((0, 1, 2), weight=1, uniform="panels")
        left = ttk.LabelFrame(columns, text="Arduino script", padding=12)
        right = ttk.LabelFrame(columns, text="Event parsing script / UDP", padding=12)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 8))
        self.parameter_editor = ParameterEditor(columns, self.reload_parameters)
        self.parameter_editor.frame.grid(row=0, column=1, sticky="nsew", padx=4)
        right.grid(row=0, column=2, sticky="nsew", padx=(8, 0))
        for frame in (left, right):
            frame.columnconfigure(1, weight=1)

        self._entry(left, "Sketch", "sketch", 0)
        button = ttk.Button(left, text="Select Arduino file…", command=self.select_sketch)
        button.grid(row=1, column=0, columnspan=3, sticky="ew", pady=(4, 10))
        self.config_widgets.append(button)
        self.port_box = self._entry(left, "Port", "port", 2, values=())
        self.refresh_button = ttk.Button(left, text="Refresh", command=lambda: self.request("refresh"))
        self.refresh_button.grid(row=2, column=2, padx=(6, 0))
        self._entry(left, "Baud", "baud", 3, values=("9600", "19200", "57600", "115200", "230400"))
        self._entry(left, "Board ID", "board", 4, values=BOARDS)
        self.detect_button = ttk.Button(left, text="Detect", command=lambda: self.request("detect", self.settings()))
        self.detect_button.grid(row=4, column=2, padx=(6, 0))
        self._entry(left, "Arduino CLI", "cli", 5)
        cli_button = ttk.Button(left, text="Browse…", command=self.select_cli)
        cli_button.grid(row=5, column=2, padx=(6, 0))
        self.config_widgets.append(cli_button)
        actions = ttk.Frame(left)
        actions.grid(row=6, column=0, columnspan=3, sticky="ew", pady=(10, 4))
        self.upload_button = ttk.Button(actions, text="Compile & Upload", command=self.upload)
        self.upload_button.pack(side="left")
        self.connect_button = ttk.Button(actions, text="Connect", command=self.toggle_connection)
        self.connect_button.pack(side="left", padx=8)
        check = ttk.Checkbutton(left, text="Connect after upload", variable=self.vars["reconnect"])
        check.grid(row=7, column=0, columnspan=3, sticky="w")
        self.config_widgets.append(check)

        self._entry(right, "Parser", "protocol", 0, protocol=True)
        picker = ttk.Button(right, text="Select event parser…", command=self.select_protocol)
        picker.grid(row=1, column=0, columnspan=3, sticky="ew", pady=(4, 10))
        self.protocol_widgets.append(picker)
        self._entry(right, "Receiver IPv4", "host", 2, protocol=True)
        self._entry(right, "UDP port", "udp_port", 3, protocol=True)
        self.protocol_button = ttk.Button(right, text="Start protocol", command=self.toggle_protocol)
        self.protocol_button.grid(row=4, column=0, columnspan=3, sticky="ew", pady=(8, 4))
        ttk.Label(right, text="Serial monitoring works without a protocol.\nUDP sends are logged; receipt is not confirmed.\nStop cancels pending sends, not the remote device state.",
                  style="Hint.TLabel", wraplength=330, justify="left").grid(row=5, column=0, columnspan=3, sticky="w", pady=8)

        command_row = ttk.LabelFrame(outer, text="Serial input", padding=10)
        command_row.grid(row=2, column=0, sticky="ew", pady=(16, 10))
        command_row.columnconfigure(0, weight=1)
        self.command_entry = ttk.Entry(command_row, textvariable=self.command)
        self.command_entry.grid(row=0, column=0, sticky="ew")
        self.command_entry.bind("<Return>", lambda event: self.send())
        ttk.Combobox(command_row, textvariable=self.vars["line_ending"], values=tuple(LINE_ENDINGS),
                     state="readonly", width=20).grid(row=0, column=1, padx=8)
        self.send_button = ttk.Button(command_row, text="Send", command=self.send)
        self.send_button.grid(row=0, column=2)

        monitor_bar = ttk.Frame(outer)
        monitor_bar.grid(row=3, column=0, sticky="ew", pady=(0, 5))
        monitor_bar.columnconfigure(0, weight=1)
        ttk.Label(monitor_bar, text="Serial monitor").grid(row=0, column=0, sticky="w")
        ttk.Checkbutton(monitor_bar, text="Autoscroll", variable=self.autoscroll).grid(row=0, column=1)
        ttk.Button(monitor_bar, text="Save log…", command=self.save_log).grid(row=0, column=2, padx=6)
        ttk.Button(monitor_bar, text="Clear", command=self.clear_log).grid(row=0, column=3)
        self.monitor = ScrolledText(outer, height=15, state="disabled", wrap="word", font=("Menlo", 11),
                                    background="#14202e", foreground="#dbe7ef", insertbackground="white")
        self.monitor.grid(row=4, column=0, sticky="nsew")
        self.monitor.tag_configure("error", foreground="#ffaaa5")
        self.monitor.tag_configure("udp", foreground="#92dac5")
        self.monitor.tag_configure("tx", foreground="#91caff")
        ttk.Label(outer, textvariable=self.status, style="Hint.TLabel").grid(row=5, column=0, sticky="w", pady=(8, 0))
        root.protocol("WM_DELETE_WINDOW", self.on_close)
        self.service.thread.start()
        self.service.log("[ready] Select a sketch and board to upload, or connect to monitor existing firmware.")
        self.request("refresh")
        self.vars["sketch"].trace_add("write", self._schedule_parameter_reload)
        root.after(50, self.poll)

    def _entry(self, frame, label, key, row, values=None, protocol=False):
        self.ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=4)
        factory = self.ttk.Entry if values is None else self.ttk.Combobox
        options = {} if values is None else {"values": values}
        widget = factory(frame, textvariable=self.vars[key], width=18, **options)
        widget.grid(row=row, column=1, sticky="ew", pady=4)
        (self.protocol_widgets if protocol else self.config_widgets).append(widget)
        return widget

    def settings(self, with_parameters=False):
        settings = Settings(**{key: var.get() for key, var in self.vars.items()})
        if with_parameters:
            if settings.sketch != self.parameter_sketch:
                if settings.sketch:
                    self._load_parameters()
                else:
                    self.parameter_editor.clear("Select an Arduino sketch to load its parameters.")
                    self.parameter_sketch = ""
            settings.parameter_overrides = self.parameter_editor.overrides()
            settings.parameter_source_hash = self.parameter_editor.hash
        return settings

    def _schedule_parameter_reload(self, *args):
        if self.parameter_reload_job:
            self.root.after_cancel(self.parameter_reload_job)
        self.parameter_reload_job = self.root.after(350, self._reload_changed_sketch)

    def _reload_changed_sketch(self):
        self.parameter_reload_job = None
        if self.vars["sketch"].get() != self.parameter_sketch:
            self.reload_parameters(quiet=True)

    def _load_parameters(self, overrides=None, expected_hash=""):
        value = self.vars["sketch"].get()
        sketch = sketch_directory(value)
        self.parameter_editor.load(sketch / (sketch.name + ".ino"), overrides, expected_hash)
        self.parameter_sketch = value

    def reload_parameters(self, quiet=False):
        from tkinter import messagebox
        try:
            self._load_parameters()
        except (OSError, ValueError) as exc:
            self.parameter_editor.clear(str(exc))
            self.parameter_sketch = None
            if not quiet:
                messagebox.showerror("Arduino parameters", str(exc), parent=self.root)

    def request(self, name, *args):
        if self.pending or self.busy or self.closing:
            return
        self.pending = True
        self.status.set("Working…")
        self.update_controls()
        self.service.submit(name, *args)

    def select_sketch(self):
        from tkinter import filedialog
        value = filedialog.askopenfilename(parent=self.root, title="Select Arduino sketch", initialdir=ROOT / "Arduino",
                                          filetypes=[("Arduino sketch", "*.ino")])
        if value:
            self.vars["sketch"].set(value)

    def select_protocol(self):
        from tkinter import filedialog
        current = self.vars["protocol"].get()
        directory = Path(current).expanduser().parent if current else PARSER_DIR
        value = filedialog.askopenfilename(parent=self.root, title="Select Python event parser",
                                          initialdir=directory, filetypes=[("Python script", "*.py")])
        if value:
            self.vars["protocol"].set(value)

    def select_cli(self):
        from tkinter import filedialog
        value = filedialog.askopenfilename(parent=self.root, title="Select Arduino CLI executable")
        if value:
            self.vars["cli"].set(value)

    def upload(self):
        from tkinter import messagebox
        try:
            settings = self.settings(with_parameters=True)
            self.request("upload", settings)
        except (OSError, ValueError) as exc:
            messagebox.showerror("Arduino parameters", str(exc), parent=self.root)

    def toggle_connection(self):
        if self.connected:
            self.request("disconnect")
        else:
            self.request("connect", self.settings())

    def toggle_protocol(self):
        if self.protocol_running:
            self.request("stop_protocol")
        else:
            self.request("start_protocol", self.settings())

    def send(self):
        if self.connected and not self.pending and not self.busy and self.command.get():
            self.request("send", self.command.get(), self.vars["line_ending"].get())

    def update_controls(self):
        blocked = self.pending or self.busy or self.closing
        self.parameter_editor.set_enabled(not blocked)
        def enable(widget, yes):
            widget.configure(state="normal" if yes else "disabled")
        for widget in self.config_widgets:
            enable(widget, not blocked and not self.connected)
        for widget in self.protocol_widgets:
            enable(widget, not blocked and not self.protocol_running)
        for widget in (self.refresh_button, self.detect_button, self.load_button):
            enable(widget, not blocked and not self.connected)
        enable(self.upload_button, not blocked)
        enable(self.connect_button, not blocked)
        enable(self.protocol_button, not blocked and self.connected)
        enable(self.send_button, not blocked and self.connected)
        enable(self.command_entry, not blocked and self.connected)
        self.connect_button.configure(text="Disconnect" if self.connected else "Connect")
        self.protocol_button.configure(text="Stop protocol" if self.protocol_running else "Start protocol")

    def poll(self):
        if self.closing:
            if not self.service.thread.is_alive():
                self.root.destroy()
            else:
                self.root.after(50, self.poll)
            return
        for _ in range(100):
            try:
                kind, value = self.service.events.get_nowait()
            except queue.Empty:
                break
            if kind == "state":
                self.connected, self.protocol_running, self.busy = value["connected"], value["protocol"], value["busy"]
                self.status.set("Compiling / uploading…" if self.busy else
                                ("Connected | Protocol running" if self.protocol_running else
                                 "Connected | Protocol stopped" if self.connected else "Disconnected"))
                self.update_controls()
            elif kind == "done":
                self.pending = False
                if value == ("send", True):
                    self.command.set("")
                self.update_controls()
            elif kind == "ports":
                self.port_box.configure(values=value)
                if not self.vars["port"].get() and len(value) == 1:
                    self.vars["port"].set(value[0])
            elif kind == "board":
                self.vars["board"].set(value)
        messages = []
        for _ in range(400):
            try:
                messages.append(self.service.logs.get_nowait())
            except queue.Empty:
                break
        if self.service.dropped:
            messages.append(f"[monitor] {self.service.dropped} display messages dropped during high traffic.")
            self.service.dropped = 0
        if messages:
            self.monitor.configure(state="normal")
            for message in messages:
                tag = "error" if "error]" in message or "failed" in message else "udp" if "[UDP" in message else "tx" if "[TX]" in message else ""
                self.monitor.insert("end", message + "\n", tag)
            count = int(self.monitor.index("end-1c").split(".")[0])
            if count > 10000:
                self.monitor.delete("1.0", f"{count - 10000 + 1}.0")
            if self.autoscroll.get():
                self.monitor.see("end")
            self.monitor.configure(state="disabled")
        self.root.after(50, self.poll)

    def save_profile(self):
        from tkinter import filedialog, messagebox
        try:
            PROFILES_DIR.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            messagebox.showerror("Save profile", str(exc), parent=self.root)
            return
        path = filedialog.asksaveasfilename(parent=self.root, initialdir=PROFILES_DIR,
                                          initialfile="behavior_profile.json", defaultextension=".json",
                                          filetypes=[("Profile", "*.json")])
        if path:
            try:
                Path(path).write_text(json.dumps(asdict(self.settings(with_parameters=True)), indent=2) + "\n", encoding="utf-8")
                self.service.log(f"[profile saved] {path}")
            except (OSError, ValueError) as exc:
                messagebox.showerror("Save profile", str(exc), parent=self.root)

    def load_profile(self):
        from tkinter import filedialog, messagebox
        path = filedialog.askopenfilename(parent=self.root, initialdir=PROFILES_DIR,
                                         filetypes=[("Profile", "*.json")])
        if path:
            try:
                data = json.loads(Path(path).read_text(encoding="utf-8"))
                defaults = asdict(Settings())
                if not isinstance(data, dict) or set(data) - set(defaults):
                    raise ValueError("Not a behavior GUI profile.")
                for key, value in data.items():
                    if type(value) is not type(defaults[key]):
                        raise ValueError(f"Invalid profile value for {key}.")
                defaults.update(data)
                if defaults["line_ending"] not in LINE_ENDINGS:
                    raise ValueError("Invalid line ending.")
                overrides = defaults["parameter_overrides"]
                if overrides:
                    sketch = sketch_directory(defaults["sketch"])
                    apply_overrides((sketch / (sketch.name + ".ino")).read_bytes(), overrides,
                                    defaults["parameter_source_hash"])
                for key in self.vars:
                    value = defaults[key]
                    self.vars[key].set(value)
                if overrides:
                    self._load_parameters(overrides, defaults["parameter_source_hash"])
                else:
                    self.reload_parameters(quiet=True)
                self.service.log(f"[profile loaded] {path}")
            except (OSError, ValueError) as exc:
                messagebox.showerror("Load profile", str(exc), parent=self.root)

    def save_log(self):
        from tkinter import filedialog, messagebox
        path = filedialog.asksaveasfilename(parent=self.root, defaultextension=".txt", filetypes=[("Text log", "*.txt")])
        if path:
            try:
                Path(path).write_text(self.monitor.get("1.0", "end-1c"), encoding="utf-8")
            except OSError as exc:
                messagebox.showerror("Save log", str(exc), parent=self.root)

    def clear_log(self):
        self.monitor.configure(state="normal")
        self.monitor.delete("1.0", "end")
        self.monitor.configure(state="disabled")

    def on_close(self):
        if self.busy or self.pending:
            self.service.log("[wait] Let the current operation finish before closing the window.")
            return
        self.closing = True
        self.status.set("Closing connections…")
        self.update_controls()
        self.service.shutdown.set()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    try:
        import tkinter as tk
    except ImportError:
        parser.exit(1, "Tkinter is required. Use a Python installation with Tk support.\n")
    root = tk.Tk()
    app = BehaviorGUI(root)
    try:
        root.mainloop()
    finally:
        app.service.shutdown.set()
        app.service.thread.join(timeout=3)


if __name__ == "__main__":
    main()
