"""Console-free desktop launcher; Conda activation is handled by conda run."""
import argparse
import ctypes
import os
from pathlib import Path
import subprocess
import tempfile
import traceback

ROOT = Path(__file__).resolve().parents[2]


def command(conda, prefix):
    return [conda, 'run', '--no-capture-output', '-p', prefix, 'python', str(ROOT / 'gui_behavior.py')]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--conda', required=True)
    parser.add_argument('--prefix', required=True)
    args = parser.parse_args()
    log_path = None
    try:
        log_dir = Path(os.environ.get('LOCALAPPDATA', tempfile.gettempdir())) / 'NeuroDAP' / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', prefix='behavior-', suffix='.log',
                                         dir=log_dir, delete=False) as log:
            log_path = Path(log.name)
            result = subprocess.run(command(args.conda, args.prefix), cwd=ROOT,
                                    stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                                    creationflags=subprocess.CREATE_NO_WINDOW)
        if result.returncode:
            raise RuntimeError(f'Behavior GUI exited with code {result.returncode}.')
        if log_path.stat().st_size == 0:
            log_path.unlink()
    except Exception as exc:
        if log_path:
            with log_path.open('a', encoding='utf-8') as log:
                log.write(traceback.format_exc())
        ctypes.windll.user32.MessageBoxW(None,
            f'{exc}\n\nDetails: {log_path or "Could not create a log file"}\n\n'
            'If you moved the repository or changed Conda, rerun gui/windows/install.py in Anaconda Prompt.',
            'NeuroDAP Behavior could not open', 0x10)


if __name__ == '__main__':
    main()
