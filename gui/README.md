# Behavior GUI

Run `gui_behavior.py` from the repository root. Its layout follows the supplied
sketch, with an additional middle column: Arduino controls at upper left,
editable Arduino parameters in the middle, event-parser/UDP controls at upper
right, then the serial input and a full-width monitor below.

## File organization

The root `gui_behavior.py` remains the launch command. Supporting files live here:

```text
gui_behavior.py          # launcher at repository root
gui/
  behavior.py           # window, serial/UDP service, and upload workflow
  parameters.py         # parameter editor and temporary sketch preparation
  manual_controls.py    # button labels and versioned serial commands
  profiles/             # default folder for saved profiles and parameter edits
  tests/                # GUI-related tests
  requirements.txt      # Python dependencies
  README.md             # this guide
```

Arduino sketches remain in `Arduino/`. The default external parser still points
to `BrainClamp/scripts/send_event_RPE.py` beside this repository.

## Setup

Use Python 3.10+ with Tkinter (included in the standard Windows/macOS Python
installers; on Linux it may need the `python3-tk` package).

```sh
python -m pip install -r gui/requirements.txt
python gui_behavior.py
```

For uploads, install [Arduino CLI](https://docs.arduino.cc/arduino-cli/installation/)
and make it available on PATH, or use the GUI's Arduino CLI Browse button. Install
the board core and any libraries required by the sketch. For an Uno, Mega, or
classic Nano, for example:

```sh
arduino-cli core update-index
arduino-cli core install arduino:avr
```

Other boards require their corresponding core. Installing Arduino IDE alone does
not guarantee that `arduino-cli` is on PATH. Compile/upload output appears in the
monitor, including missing-core or missing-library errors. Board discovery is
supported by [Arduino CLI](https://docs.arduino.cc/arduino-cli/getting-started/);
unknown boards can still be selected manually.

## Workflow

1. **Select Arduino file**: choose a `.ino`. The main file must have the same name
   as its parent folder, as required by Arduino. The entire sketch folder is
   compiled, including accompanying tabs and headers. Selection alone does not
   upload anything. One existing repository folder,
   `Arduino/Shun_OptoPair_EPLHbPair_random`, does not match its main `.ino` name;
   arrange a correctly named copy before uploading that sketch.
2. **Port / Board ID**: refresh ports, select the connected Arduino, then use
   Detect or select/type a board ID. For a Mega 2560 use
   `arduino:avr:mega:cpu=atmega2560`. The defaults are COM4, 115200 baud, and this
   Mega 2560 board ID. You can change them before connecting.
3. **Arduino parameters / Compile & Upload**: selecting a sketch fills the middle
   column from its marked user settings. Edit `UnitRewardSize` in the pinned top
   row for daily valve calibration; scroll through the other settings below it.
   Click Compile & Upload to validate your edits and copy the entire sketch to a
   temporary folder. The app applies edits to that copy, stops the protocol,
   cancels pending UDP work, closes serial, compiles, and uploads only on compile
   success. Original sketch files are never rewritten. The monitor shows the
   parameter overrides and upload progress. Controls are disabled during this work.
   With **Connect after upload** checked, it attempts to reopen the selected
   port. Otherwise click Connect yourself. Upload failures leave the app
   disconnected. No protocol is automatically restarted.
4. **Serial input**: enter a command and press Enter or Send. Select LF, CR, CRLF,
   or no line ending to match the sketch. Commands are UTF-8; whitespace is
   preserved. Successful sends clear the textbox; failed sends retain the text.
   Set the baud rate to the value in the sketch's `Serial.begin(...)`.
5. **Optional event protocol**: select the Python parser, set the receiver IPv4
   address and UDP port, then click Start protocol **before opening
   `brainclamp_gui.py`**. It handles new serial lines
   while you continue using the command box. Stop protocol keeps serial connected.
   The default parser is `BrainClamp/scripts/send_event_RPE.py`, with BrainClamp
   next to the NeuroDAP repository. For example, on your Windows rig this resolves
   to `C:/Shun-local/BrainClamp/scripts/send_event_RPE.py`. The parser picker opens
   the selected parser's folder. Nothing runs until you explicitly start it.

Close any other serial monitor before connecting. Opening a serial connection can
reset an Arduino. Allow startup to finish before sending commands or starting the
experiment. Some boards change port names after upload; refresh and reconnect if
automatic reconnection fails. This app does not send a task-start command on its
own: use Start task or the command understood by your firmware.

**Stop protocol stops local forwarding and cancels delayed/repeated packets. It
does not undo commands already received by the other computer or stop the Arduino
experiment.** Use the device's own stop controls when you need to stop it.

## Manual controls and stopping

The shorter configuration panels leave room for grouped manual buttons above the
free-text serial input. Upload the updated **Shun_DAClamp_Reward** or
**Shun_DAClamp_Random** sketch to enable them. Each sketch includes a local
`BehaviorManual.h`; keep that file beside its `.ino` when copying the sketch.
The GUI waits for `BEHAVIOR_CONTROLS 1` from the Arduino before enabling buttons,
so an older sketch cannot interpret a button using a conflicting digit mapping.
It queries once after connecting; **Check controls** sends another query if needed.
Other sketches still support the free-text serial monitor.

The combined **Lasers** group contains blue and red open/close toggles plus one
pattern button per color. The **Session** group has Water calibration and one
large task toggle: green **Start task**, red **End task**. The protocol toggle
uses the same green Start / red Stop colors. Task and laser labels follow Arduino
acknowledgments, including commands sent through the textbox; rejected commands
do not change their state. Updated firmware also reports its current task and
manual laser state on startup and when queried by **Check controls**.
During water calibration the task toggle also shows red **End task**, so you can
cancel calibration immediately; it returns to Start when calibration finishes.

| Button | Serial byte | Behavior |
| --- | --- | --- |
| Small reward | `w` | Water for `SmallRewardSize` ms |
| Large reward | `W` | Water for `BigRewardSize` ms |
| Punishment (no tone) | `p` | Airpuff for `SmallPunishSize` ms; does not start a tone |
| Tone | `t` | `LeftCueFreq` for `ShortToneDuration` |
| Open / close blue | `b` / `B` | Hold blue shutter open / close it and cancel its pattern |
| Open / close red | `r` / `R` | Hold red shutter open / close it and cancel its pattern |
| Blue / red pattern | `f` / `F` | Run the corresponding pattern on the Arduino |
| Water calibration | `c` | Deliver `CalibrationRepeats` unit rewards |
| Start task | `s` | Start the task at its ITI with existing counts retained |
| End task | `x` | Stop task and outputs, preserving recorded counts |

**End task does not reset counts or reboot the board.** It closes both water
valves and the airpuff, closes both shutters, stops tone, and cancels patterns and
calibration. Trial, outcome, manual reward/punishment, cue, and performance counts
are retained; the completed-calibration delivery count is also retained. Repeated
End requests are harmless. A later Start begins another ITI with the retained
session counts rather than resuming an interrupted trial. Uploading firmware or
resetting the physical board still reinitializes RAM counters, as before.

The GUI's End button also stops local Python event forwarding and cancels queued
UDP sends. It does not undo commands already received by a remote computer.
Start task does not automatically start the Python protocol; enable that separately
when needed. Each manual request is logged, followed by the Arduino's `ACK`,
`DONE`, or `ERR` response. `ACK` confirms command acceptance, not a physical
measurement of delivered reward or light.

Manual stimuli and calibration require the task to be stopped, preventing them
from competing with scheduled task outputs. Close-shutter and End commands remain
available during a task or calibration. Starting a task closes any manual outputs;
an active calibration must finish or be stopped first.

Both sketches now expose `BluePatternPulseMs`, `BluePatternPeriodMs`,
`BluePatternPulseCount` and equivalent red parameters in the middle panel. Initial
defaults are **25 pulses, 5 ms wide, 50 ms period (20 Hz)** per color. The period
must be at least the pulse width, and width/count must be positive. The Arduino
generates patterns using elapsed-time scheduling; the GUI only starts them.
Check these parameters for your experiment before uploading.

Calibration defaults to **200 deliveries** with **100 ms between deliveries**.
It uses the uploaded `UnitRewardSize` and runs without blocking serial processing,
so End can close the valve and cancel remaining deliveries immediately on receipt.
The Random task's tone also uses nonblocking timing so it can be interrupted.

Legacy digits remain available in the textbox: `1` is small reward in Reward and
large reward in Random; `2` is tone in Reward and punishment in Random; `3`/`4`
open shutters, `5`/`6` run patterns, `7` calibrates, and `8`/`9` start/stop.
The same idle-only restriction applies to legacy manual-stimulus commands.

## Editing Arduino parameters

The middle panel reads the main `.ino` file's **User settings** section, ending
at **Params Initializtion** (or **Params Initialization**). These markers are
already present in the repository's behavioral-task sketches. `UnitRewardSize`
is pinned above the scrolling area and is entered as a positive whole number of
milliseconds. This is the valve-open duration for one calibrated reward unit;
actual volume depends on your rig's calibration.

Other supported settings include scalar numbers, `true`/`false` booleans, simple
arithmetic expressions, and fixed-size one-dimensional arrays. Original names,
section headings, and inline comments are displayed. For example:

- Changing `UnitRewardSize` from `24` to `30` preserves
  `SmallRewardSize = 2 * UnitRewardSize`, which then compiles to 60 ms.
- You can edit `BigRewardSize` to `8 * UnitRewardSize` to change its multiplier.
- Edit a probability range using braces, such as `{1, 40}`.

Values are source initializers, not live readbacks from the Arduino. Some sketches
override initial values later in `setup()`; the GUI preserves that firmware logic.
Editing while connected is allowed, but takes effect only after another upload.
**Reload / reset edits** discards the current edits and rereads the
source. Selecting another sketch loads its own parameters without carrying edits
over from the previous sketch.
The adjacent **Compile & Upload** button performs exactly the same operation as
the button in the Arduino script panel, including applying the edited parameters.

For another sketch, mark the intended configuration region explicitly:

```cpp
// GUI parameters begin
unsigned long UnitRewardSize = 24;
unsigned long SmallRewardSize = 2 * UnitRewardSize;
boolean ENL = true;
int PairProbRange[2] = {1, 40};
// GUI parameters end
```

Put one initialized declaration per line. The editor intentionally skips runtime
state outside the markers, comments, function-local variables, conditional
preprocessor blocks, duplicate names, `#define`s, and complex declarations or
initializers (such as function calls). A sketch without a supported marked region
can still be uploaded normally; its parameter panel explains that no settings
were found. The pinned field reads **Not defined** if `UnitRewardSize` is absent.

Invalid edits are rejected before compiling. Arduino CLI still performs the final
C++ type/compile checks. The GUI records a fingerprint of the main source when
loading parameters. If that source changes on disk before you upload edited
values, it asks you to reload and reapply them rather than silently overwriting
new source values. Files within the sketch folder (including headers, extra tabs,
and `src/`) are copied for compilation; keep local sketch dependencies inside
that folder rather than referencing files outside it with `../` includes.

## External protocol scripts

Protocol scripts remain in your external `BrainClamp/scripts` folder; no copy is
bundled with the GUI. The default is `send_event_RPE.py`, and Select event parser
lets you choose another script from that folder or elsewhere.

The loader imports
its event functions without calling `main()` or constructing its GUI/serial
connection. The GUI supplies the UDP sender and overrides its `delay_send()`
helper with cancellable scheduling. Its hard-coded COM port, UDP address, baud,
`FORWARD_EVENTS`, and `LOG_JSONL` settings are not used. Each Start loads a fresh
module, resetting protocol globals and applying source edits.

UDP format remains `MSG <sequence> <command>\n`, with three copies per command
and a nominal 3 ms gap. Sequence IDs continue across protocol restarts in this
GUI session and begin at a timestamp-derived value to avoid immediate collisions
with a previous app session. The BrainClamp receiver must deduplicate IDs.
The monitor reports **sent**, not acknowledged delivery. Scheduling uses the PC
clock and worker loop; it is not a guarantee of millisecond-accurate delivery.

## Other event scripts

Select a trusted local Python module defining these two functions:

```python
def parse_arduino_line(line):
    # Return any event object, or None to ignore the line.
    return {"type": "cue"} if line.startswith("Cue") else None

def apply_session_logic(event, udp):
    if event["type"] == "cue":
        udp.send("CMD PID_ON")
        udp.call_later("CMD PID_OFF", 0.5)
```

Keep GUI startup and serial-opening code behind `if __name__ == "__main__":`.
Module-level code executes on import; these scripts are ordinary trusted Python,
not a sandbox. Do not open the same serial port from the protocol. Event callbacks
should return promptly. Use `udp.call_later(command, seconds)` for delayed work
instead of sleeping, starting threads, or creating independent sockets. The
legacy `delay_send(udp, command, seconds)` helper is supported as described above.
Protocol errors stop forwarding and cancel pending sends while the raw monitor
continues. Slow callbacks can delay serial processing and stop actions.

## Profiles and logs

Save/load profiles to remember sketch, board, port, parser, baud, line ending,
CLI location, UDP settings, and Arduino parameter overrides (including daily
`UnitRewardSize` calibration). Profiles store only changed initializers and a
source fingerprint; they do not contain a copy of the sketch. Older profiles
without parameter overrides still work. If a profile's source fingerprint no
longer matches the selected sketch, its parameter overrides are rejected.
**Save profile** opens a file dialog in `gui/profiles/`: choose a filename for a
`.json` file or browse to another folder. Parameter edits are stored in that same
file, not in a separate parameter file. There is no automatic saving. The monitor
prints the saved file's full path. **Load
profile** opens a previously saved `.json` file and restores its settings into
the controls; its file dialog also starts in `gui/profiles/`. Existing profiles
saved elsewhere can still be selected. Profiles are not loaded automatically at startup; the app starts
with the defaults above. Loading a profile never connects or sends data.
Profiles contain local file paths; update them when moving between computers.

The monitor keeps the most recent 10,000 lines, showing raw RX, TX, parsed event
types, UDP sends, and upload output. Save log exports the currently displayed
text. Under heavy traffic, display messages may be dropped with an explicit
notice; this does not skip event parsing. This is a monitor, not a complete
experiment-data recorder, and does not duplicate the original JSONL EventStore.

## Verification

```sh
python -m unittest discover -s gui/tests -p 'test_behavior*.py' -v
```

Tests use simulated serial, UDP, and CLI components; no Arduino or remote
experiment is operated. When a C++ compiler is available, both complete DA-clamp
sketches are compiled and executed against mock Arduino I/O to check stop/count
preservation, patterns, rewards, calibration interruption and legacy commands.
This does not replace an Arduino-core build or validation on the target rig.
