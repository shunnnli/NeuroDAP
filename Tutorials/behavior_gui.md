# Behavior GUI

Run `behavior_gui.py` from the repository root. Its layout follows the supplied
sketch: Arduino controls at upper left, event-parser/UDP controls at upper right,
then the serial input and a full-width monitor below.

## Setup

Use Python 3.10+ with Tkinter (included in the standard Windows/macOS Python
installers; on Linux it may need the `python3-tk` package).

```sh
python -m pip install -r requirements-behavior.txt
python behavior_gui.py
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
   `arduino:avr:mega:cpu=atmega2560`. Confirm the detected selection. Board ID is
   left blank initially so the app does not assume your hardware.
3. **Compile & Upload**: the app stops the protocol, cancels pending UDP work,
   closes serial, compiles to a temporary directory, and uploads only on compile
   success. The monitor shows progress. Controls are disabled during this work.
   With **Connect after upload** checked, it attempts to reopen the selected
   port. Otherwise click Connect yourself. Upload failures leave the app
   disconnected. No protocol is automatically restarted.
4. **Serial input**: enter a command and press Enter or Send. Select LF, CR, CRLF,
   or no line ending to match the sketch. Commands are UTF-8; whitespace is
   preserved. Successful sends clear the textbox; failed sends retain the text.
   Set the baud rate to the value in the sketch's `Serial.begin(...)`.
5. **Optional event protocol**: select the Python parser, set the receiver IPv4
   address and UDP port, then click Start protocol. It handles new serial lines
   while you continue using the command box. Stop protocol keeps serial connected.
   The bundled `behavior_protocols/rpe.py` is selected initially, but does not run
   until you explicitly start it.

Close any other serial monitor before connecting. Opening a serial connection can
reset an Arduino. Allow startup to finish before sending commands or starting the
experiment. Some boards change port names after upload; refresh and reconnect if
automatic reconnection fails. This app does not send a task-start command on its
own: use the command understood by your firmware.

**Stop protocol stops local forwarding and cancels delayed/repeated packets. It
does not undo commands already received by the other computer or stop the Arduino
experiment.** Use the device's own stop controls when you need to stop it.

## RPE compatibility

The bundled protocol preserves the active logic in the supplied
`BrainClamp/scripts/send_event_RPE.py`:

| Arduino event | UDP behavior |
| --- | --- |
| Cue / `TRIAL_ON` | `PID_ON 200`, then after 0.4 s `SET_TARGET 0.8000 0.25 natural` |
| Big or small reward | `SET_TARGET 0.8000`, then after 0.5 s `PID_OFF 500` |
| Omission | `SET_TARGET -0.1000 0.25 natural`, then after 0.5 s `PID_OFF 500` |
| Other outcome or timeout | After 0.5 s `PID_OFF 500` |

The parser also recognizes ITI, licks, task/test banners, manual rewards, and
OperantClamp DA clamp/unclamp messages. Those events do not send commands unless
the protocol logic explicitly handles them. A different Arduino sketch may print
different messages and need a different parser.

You may select the **original `send_event_RPE.py` directly**. The loader imports
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
CLI location, and UDP settings. Loading a profile never connects or sends data.
Profiles contain local file paths; update them when moving between computers.

The monitor keeps the most recent 10,000 lines, showing raw RX, TX, parsed event
types, UDP sends, and upload output. Save log exports the currently displayed
text. Under heavy traffic, display messages may be dropped with an explicit
notice; this does not skip event parsing. This is a monitor, not a complete
experiment-data recorder, and does not duplicate the original JSONL EventStore.

## Verification

```sh
python -m unittest discover -s tests -p 'test_behavior_gui.py' -v
```

Tests use simulated serial, UDP, and CLI components; no Arduino or remote
experiment is operated. Real board upload and receiver behavior still require
verification on the target rig.
