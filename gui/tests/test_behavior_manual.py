"""GUI command tests and execution of both sketches against mocked Arduino I/O."""
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import Mock

from gui.behavior import BehaviorService, ROOT
from gui.manual_controls import CAPABILITY_LINE, COMMANDS


class ManualCommandTests(unittest.TestCase):
    def test_controls_require_firmware_identification(self):
        service = BehaviorService()
        service.ser = Mock()
        with self.assertRaisesRegex(RuntimeError, 'firmware'):
            service.do_manual('small_reward')
        service.ser.write.assert_not_called()
        packet = (CAPABILITY_LINE + '\n').encode()
        service.ser.read.return_value = packet
        service.ser.in_waiting = len(packet)
        service._read_serial()
        self.assertTrue(service.manual_ready)

    def test_commands_send_one_byte_and_end_cancels_udp_without_reset_command(self):
        service = BehaviorService()
        service.ser = Mock()
        service.ser.write.return_value = 1
        service.manual_ready = True
        for action, command in COMMANDS.items():
            service.do_manual(action)
            service.ser.write.assert_called_with(command.encode())
        service.protocol = Mock()
        protocol = service.protocol
        service.udp = Mock()
        udp = service.udp
        service.do_manual('end_task')
        udp.close.assert_called_once()
        protocol.close.assert_called_once()
        service.ser.write.assert_called_with(b'x')
        service.do_disconnect()
        self.assertFalse(service.manual_ready)


@unittest.skipUnless(shutil.which('clang++') or shutil.which('g++'), 'C++ compiler unavailable')
class FirmwareTests(unittest.TestCase):
    def test_both_sketches_keep_counts_on_stop_and_execute_manual_actions(self):
        compiler = shutil.which('clang++') or shutil.which('g++')
        flags = []
        if sys.platform == 'darwin':
            sdk = subprocess.check_output(['xcrun', '--show-sdk-path'], text=True).strip()
            flags = ['-isysroot', sdk, '-isystem', str(Path(sdk) / 'usr/include/c++/v1')]
        harness = Path(__file__).parent / 'firmware'
        headers = []
        for name in ('Shun_DAClamp_Reward', 'Shun_DAClamp_Random'):
            with self.subTest(sketch=name), tempfile.TemporaryDirectory() as temp:
                sketch = ROOT / 'Arduino' / name / (name + '.ino')
                headers.append((sketch.parent / 'BehaviorManual.h').read_bytes())
                # Arduino normally generates these prototypes before C++ compilation.
                declarations = re.findall(r'\b(void|bool)\s+(\w+)\(([^)]*)\)\s*\{', sketch.read_text())
                source = '#include "Arduino.h"\n'
                source += '\n'.join(f'{kind} {func}({args});' for kind, func, args in declarations)
                source += f'\n#include "{sketch.as_posix()}"\n'
                source += f'#include "{(harness / "manual_controls_checks.cpp").resolve().as_posix()}"\n'
                path = Path(temp) / 'check.cpp'
                path.write_text(source)
                executable = Path(temp) / 'check'
                built = subprocess.run([compiler, '-std=c++11', *flags, '-I', str(harness), str(path), '-o', str(executable)],
                                       capture_output=True, text=True, timeout=60)
                self.assertEqual(built.returncode, 0, built.stderr)
                ran = subprocess.run([str(executable)], capture_output=True, text=True, timeout=10)
                self.assertEqual(ran.returncode, 0, ran.stderr)
        self.assertEqual(headers[0], headers[1], 'The two standalone manual-control headers must match')
