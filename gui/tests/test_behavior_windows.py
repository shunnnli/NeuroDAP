"""Portable checks for environment selection and the Windows launch command."""
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from gui.windows import install, launch


class WindowsSetupTests(unittest.TestCase):
    def test_reuses_working_brainclamp_without_creating_environment(self):
        with patch.object(install, 'run_json', return_value={'envs': ['/conda/envs/brainclamp']}), \
             patch.object(install, 'runtime_works', return_value=True), \
             patch.object(install.subprocess, 'run') as run:
            self.assertEqual(install.select_environment('conda'), '/conda/envs/brainclamp')
            run.assert_not_called()

    def test_failed_brainclamp_check_falls_back_without_modifying_it(self):
        with patch.object(install, 'run_json', return_value={'envs': ['/envs/brainclamp', '/envs/neurodap']}), \
             patch.object(install, 'runtime_works', side_effect=[False, True]), \
             patch.object(install.subprocess, 'run') as run:
            self.assertEqual(install.select_environment('conda'), '/envs/neurodap')
            run.assert_not_called()

    def test_creates_neurodap_only_when_needed(self):
        with patch.object(install, 'run_json', side_effect=[{'envs': []}, {'envs': ['/envs/neurodap']}]), \
             patch.object(install, 'runtime_works', return_value=True), \
             patch.object(install.subprocess, 'run') as run:
            self.assertEqual(install.select_environment('conda'), '/envs/neurodap')
            run.assert_called_once_with(['conda', 'create', '-y', '-n', 'neurodap',
                                         'python=3.11', 'tk', 'pip'], check=True)

    def test_explicit_missing_environment_is_not_silently_replaced(self):
        with patch.object(install, 'run_json', return_value={'envs': []}), \
             patch.object(install.subprocess, 'run') as run:
            with self.assertRaises(RuntimeError):
                install.select_environment('conda', 'brainclamp')
            run.assert_not_called()

    def test_launch_keeps_paths_as_arguments_and_runs_current_repo(self):
        command = launch.command(r'C:\User Files\Conda\Scripts\conda.exe', r'C:\User Files\envs\brainclamp')
        self.assertEqual(command[4], r'C:\User Files\envs\brainclamp')
        self.assertEqual(command[-1], str(launch.ROOT / 'gui_behavior.py'))
        self.assertEqual(command[1:4], ['run', '--no-capture-output', '-p'])

    def test_shortcut_paths_are_data_not_powershell_source(self):
        with patch.object(install, 'run_json', return_value={'root_prefix': '/Conda With Spaces'}), \
             patch.object(install, 'shortcut_icon', return_value=Path('/icons/behavior-new.ico')), \
             patch.object(Path, 'is_file', return_value=True), \
             patch.object(install.subprocess, 'run') as run:
            install.create_shortcut('/Conda With Spaces/conda.exe', '/envs/brainclamp')
            args, kwargs = run.call_args
            self.assertEqual(args[0][-1], install.SHORTCUT_SCRIPT)
            self.assertIn('NEURODAP_SHORTCUT_DATA', kwargs['env'])
            self.assertTrue(kwargs['check'])

    def test_icon_revision_changes_path_only_when_artwork_changes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / 'gui').mkdir()
            source = root / 'gui/icon-windows.ico'
            source.write_bytes(b'first icon')
            with patch.object(install, 'ROOT', root), \
                 patch.dict(install.os.environ, {'LOCALAPPDATA': str(root / 'Local App Data')}):
                first = install.shortcut_icon()
                self.assertEqual(first.read_bytes(), source.read_bytes())
                self.assertEqual(install.shortcut_icon(), first)
                source.write_bytes(b'updated icon')
                second = install.shortcut_icon()
                self.assertNotEqual(first, second)
                self.assertEqual(second.read_bytes(), source.read_bytes())
                self.assertTrue(first.exists())  # other shortcuts may still refer to it


if __name__ == '__main__':
    unittest.main()
