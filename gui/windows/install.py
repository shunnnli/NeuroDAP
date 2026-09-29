r"""Run in Anaconda Prompt: python gui\windows\install.py."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
CHECK_RUNTIME = "import sys, tkinter; assert sys.version_info >= (3, 10); r=tkinter.Tk(); r.withdraw(); r.update(); r.destroy()"
SHORTCUT_SCRIPT = r'''
$ErrorActionPreference = 'Stop'
$data = $env:NEURODAP_SHORTCUT_DATA | ConvertFrom-Json
$shell = New-Object -ComObject WScript.Shell
$path = Join-Path ($shell.SpecialFolders.Item('Desktop')) 'NeuroDAP Behavior.lnk'
$link = $shell.CreateShortcut($path)
$link.TargetPath = $data.target
$link.Arguments = $data.arguments
$link.WorkingDirectory = $data.directory
$link.IconLocation = $data.target + ',0'
$link.Description = 'NeuroDAP behavior control (live repository version)'
$link.Save()
Write-Output $path
'''


def run_json(args):
    return json.loads(subprocess.check_output(args, text=True, encoding='utf-8'))


def find_conda():
    candidates = [os.environ.get('CONDA_EXE'), shutil.which('conda.exe')]
    # Anaconda Prompt normally supplies CONDA_EXE, even when conda is a shell wrapper.
    for candidate in candidates:
        if candidate and Path(candidate).is_file() and Path(candidate).suffix.lower() == '.exe':
            return str(Path(candidate).resolve())
    raise RuntimeError('Open Anaconda Prompt (or Miniforge Prompt) and run this installer there.')


def named_environment(envs, name):
    return next((p for p in envs if Path(p).name.lower() == name.lower()), None)


def runtime_works(conda, prefix):
    result = subprocess.run([conda, 'run', '-p', prefix, 'python', '-c', CHECK_RUNTIME],
                            capture_output=True, text=True)
    return result.returncode == 0


def select_environment(conda, requested=None):
    envs = run_json([conda, 'env', 'list', '--json'])['envs']
    if requested:
        prefix = named_environment(envs, requested)
        if not prefix or not runtime_works(conda, prefix):
            raise RuntimeError(f'Environment {requested!r} must exist and provide Python 3.10+ and working Tkinter.')
        return prefix
    brainclamp = named_environment(envs, 'brainclamp')
    if brainclamp and runtime_works(conda, brainclamp):
        return brainclamp
    if brainclamp:
        print('BrainClamp did not pass the Python/Tkinter check; using neurodap instead.')
    prefix = named_environment(envs, 'neurodap')
    if not prefix:
        subprocess.run([conda, 'create', '-y', '-n', 'neurodap', 'python=3.11', 'tk', 'pip'], check=True)
        prefix = named_environment(run_json([conda, 'env', 'list', '--json'])['envs'], 'neurodap')
    if not prefix or not runtime_works(conda, prefix):
        raise RuntimeError('The neurodap environment needs Python 3.10+ and working Tkinter. Repair it and rerun setup.')
    return prefix


def create_shortcut(conda, prefix):
    base = Path(run_json([conda, 'info', '--json'])['root_prefix'])
    pythonw = base / 'pythonw.exe'
    if not pythonw.is_file():
        raise RuntimeError(f'Cannot find the Conda launcher: {pythonw}')
    data = {'target': str(pythonw), 'directory': str(ROOT),
            'arguments': subprocess.list2cmdline([str(ROOT / 'gui/windows/launch.py'),
                                                '--conda', conda, '--prefix', prefix])}
    env = os.environ.copy()
    env['NEURODAP_SHORTCUT_DATA'] = json.dumps(data)
    subprocess.run(['powershell.exe', '-NoProfile', '-NonInteractive', '-Command', SHORTCUT_SCRIPT],
                   env=env, check=True)


def main():
    parser = argparse.ArgumentParser(description='Install the NeuroDAP desktop shortcut on Windows.')
    parser.add_argument('--env', help='Use this existing Conda environment instead of automatic selection.')
    args = parser.parse_args()
    if sys.platform != 'win32':
        parser.exit(1, 'This installer is for Windows.\n')
    try:
        conda = find_conda()
        prefix = select_environment(conda, args.env)
        print(f'Using environment: {prefix}', flush=True)
        subprocess.run([conda, 'run', '--no-capture-output', '-p', prefix, 'python', '-m', 'pip',
                        'install', '-r', str(ROOT / 'gui/requirements.txt')], check=True)
        subprocess.run([conda, 'run', '-p', prefix, 'python', '-c',
                        'import serial; from gui.behavior import BehaviorGUI'], cwd=ROOT, check=True)
        create_shortcut(conda, prefix)
        print('Ready: double-click NeuroDAP Behavior on your desktop. No manual activation needed.')
        if not (ROOT.parent / 'BrainClamp/scripts/send_event_RPE.py').is_file():
            print('Parser not found beside NeuroDAP. Select your BrainClamp parser in the GUI.')
        if not shutil.which('arduino-cli'):
            print('Arduino CLI is not on PATH. See gui/windows/README.md for the one-time upload setup.')
    except (RuntimeError, OSError, subprocess.SubprocessError, ValueError, KeyError) as exc:
        parser.exit(1, f'Setup failed: {exc}\nFix the issue and rerun the same installer.\n')


if __name__ == '__main__':
    main()
