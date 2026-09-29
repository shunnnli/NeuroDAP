# Windows desktop setup

Recommended: reuse the existing **brainclamp** Conda environment. The installer
checks Python 3.10+ and Tkinter, then installs the small GUI requirements list
(currently pySerial). It does not reinstall BrainClamp or its scientific packages.
If brainclamp is missing or fails the runtime check, it uses or creates a lightweight
**neurodap** environment with Python 3.11, Tkinter, pip, and the GUI requirements.
An existing but broken neurodap environment is reported rather than overwritten.

## On each computer

1. Install Conda and clone NeuroDAP. If BrainClamp is already installed, keep it.
   Put the repositories beside each other for automatic parser discovery:

   ```text
   C:\Shun-local\NeuroDAP\
   C:\Shun-local\BrainClamp\
   ```

   Other parent folders work too. If BrainClamp is elsewhere, use **Select event
   parser** in the GUI and save a profile; load that profile on future launches.

2. Open **Anaconda Prompt** or **Miniforge Prompt** and run:

   ```bat
   cd /d C:\Shun-local\NeuroDAP
   python gui\windows\install.py
   ```

   Change the first path to your actual repository folder. No manual environment
   creation or activation is needed. Internet access is needed to create an
   environment or install missing packages. If Conda asks you to resolve channel
   terms or configuration, do that in the prompt and rerun setup.

3. Double-click **NeuroDAP Behavior** on your desktop. The shortcut uses the supplied
   NeuroDAP icon, runs the current repository code, and handles Conda activation invisibly.
   It does not start BrainClamp, upload firmware, or start a protocol automatically.

To choose an existing environment explicitly, append `--env brainclamp` or
`--env neurodap` to the installer command. An explicitly selected environment must
already exist and pass the runtime check. Environment sharing does not prevent
BrainClamp and this GUI from running at the same time; only one application may
own the Arduino serial port.

## Arduino uploads on a new computer

Skip this if Arduino CLI and the AVR board core already work on the computer.
Otherwise, install [Arduino CLI](https://docs.arduino.cc/arduino-cli/installation/).
With Windows Package Manager available, run in a terminal:

```bat
winget install --id ArduinoSA.CLI --exact --source winget
```

Close and reopen the terminal after installation so it sees the updated PATH,
then install the core used by the Mega, Uno, and classic Nano:

```bat
arduino-cli core update-index
arduino-cli core install arduino:avr
```

These steps install tools only; they do not upload anything to the Arduino.
If CLI is installed outside PATH, select its executable using the GUI's Browse
button and save a profile. Other boards/sketches may need different cores or
libraries. USB drivers may also be needed for the particular board.

## Updates and troubleshooting

- To apply the custom icon to a shortcut created before icon support was added,
  rerun `python gui\windows\install.py` once. The installer uses
  `gui/icon-windows.ico`, converted from `gui/icon-windows.png` with multiple sizes.
  Restart the GUI for its window icon to update as well.

- Pull NeuroDAP updates normally, then close and reopen the GUI. The shortcut
  points to the repository; no app rebuild or shortcut recreation is needed.
- If requirements change, rerun the installer. If you move the repository,
  relocate Conda, or remove its environment, rerun it to update the shortcut.
- The shortcut is stored on the Windows Desktop (including redirected/OneDrive
  desktops). It can be moved elsewhere after creation. Rerunning setup recreates
  the Desktop copy with the same name.
- Startup failures show a message and save a log under
  `%LOCALAPPDATA%\NeuroDAP\logs`. Failed startup and nonempty output logs are kept;
  empty successful logs are removed.
- The installer does not download BrainClamp or provision its hardware drivers.
  Its supplied RPE parser uses pySerial and the Python standard library; other
  parsers may need additional packages in the chosen environment.
- Remove the shortcut to remove the desktop entry. No registry or startup entries
  are created, and no Conda environment is deleted.

The launcher uses [conda run](https://docs.conda.io/projects/conda/en/stable/commands/run.html)
so environment activation hooks and Windows DLL paths are applied automatically.
