"""RPE event protocol adapted from scripts/send_event_RPE.py in BrainClamp.

Loaded by behavior_gui.py. No serial connection or GUI is created here.
Edit apply_session_logic() to change the experimental protocol.
"""
from __future__ import annotations

import re
import time
from dataclasses import dataclass
from typing import Any, Dict, Optional


def send_protocol(udp, cmd, rep=None, rep_max=50):
    return udp.send(cmd)


def delay_send(udp, cmd, delay_s):
    if udp is not None:
        udp.call_later(cmd, delay_s)


@dataclass
class Event:
    type: str                 # e.g. "TRIAL_ON", "ITI", "OUTCOME", "TIMEOUT", "TEST_START"
    wall_time_s: float        # time.time() on the event PC when we parsed it
    trial: Optional[int]      # TrialNum if available
    arduino_time_s: Optional[float]  # time printed by Arduino (seconds), if present
    raw: str                  # original serial line
    meta: Dict[str, Any]


TRIAL_PREFIX = re.compile(r"^Trial:\s*(\d+)\s+")


FLOAT_AFTER = lambda key: re.compile(rf"{re.escape(key)}\s*:\s*([0-9]*\.?[0-9]+)")


RE_ITI = re.compile(r"Trial:\s*(\d+)\s+ITI\s*=\s*([0-9]*\.?[0-9]+)\s+Time:\s*([0-9]*\.?[0-9]+)")


RE_CUE = re.compile(
    r"Trial:\s*(\d+)\s+Cue start\s*\((?P<label>[^#\)]+)"
    r"(?:\s*#(?P<idx>\d+))?\)\s+Time:\s*(?P<t>[0-9]*\.?[0-9]+)"
)


RE_OUTCOME = re.compile(
    r"Trial:\s*(\d+)\s+Reward block:\s*(?P<label>.+?)\s+Time:\s*(?P<t>[0-9]*\.?[0-9]+)"
)


RE_TIMEOUT = re.compile(
    r"Trial:\s*(\d+)\s+Time out\s+Time:\s*(?P<t>[0-9]*\.?[0-9]+)"
)


RE_TEST_START = re.compile(r"^\*{5,}New test started\*{5,}$")


RE_TEST_END = re.compile(r"^\*{5,}Stim test finished\*{5,}$")


RE_LICK = re.compile(
    r"Trial:\s*(\d+)\s+Lick Detected\s+Time:\s*(?P<t>[0-9]*\.?[0-9]+)"
)


RE_TASK_START = re.compile(r"TASK STARTED AT\s+([0-9]*\.?[0-9]+)")


RE_TASK_END = re.compile(r"TASK ENDED AT\s+([0-9]*\.?[0-9]+)")


RE_MANUAL_REWARD = re.compile(r"Manual reward\s+(\d+)\s+(\d+)")


RE_MANUAL_PUNISH = re.compile(r"Manual punishment\s+(\d+)\s+(\d+)")


RE_SIMPLE_ITI = re.compile(r"^ITI:\s*([0-9]*\.?[0-9]+)")


RE_SIMPLE_LICK = re.compile(r"Lick Detected.*Time:\s*([0-9]*\.?[0-9]+)")


RE_DA_UNCLAMP = re.compile(
    r"DA unclamped\s*#(?P<count>\d+)\s*Time:\s*(?P<t>[0-9]*\.?[0-9]+)"
)


RE_DA_CLAMP = re.compile(
    r"DA reclamped\s*#(?P<count>\d+)\s*Time:\s*(?P<t>[0-9]*\.?[0-9]+)"
)


RE_MANUAL_DA_UNCLAMP = re.compile(
    r"Manual DA unclamp\s*#(?P<count>\d+)\s*Time:\s*(?P<t>[0-9]*\.?[0-9]+)"
)


RE_MANUAL_DA_CLAMP = re.compile(
    r"Manual DA clamp\s*#(?P<count>\d+)\s*Time:\s*(?P<t>[0-9]*\.?[0-9]+)"
)


def parse_arduino_line(line: str) -> Optional[Event]:
    """
    Parse your existing Serial prints into structured events.
    Returns None if not recognized.
    """
    s = line.strip()
    now = time.time()

    # Test banners
    if RE_TEST_START.match(s):
        return Event(type="TEST_START", wall_time_s=now, trial=None, arduino_time_s=None, raw=s, meta={})
    if RE_TEST_END.match(s):
        return Event(type="TEST_END", wall_time_s=now, trial=None, arduino_time_s=None, raw=s, meta={})

    # ITI line
    m = RE_ITI.search(s)
    if m:
        trial = int(m.group(1))
        iti_s = float(m.group(2))
        t_s = float(m.group(3))
        return Event(
            type="ITI",
            wall_time_s=now,
            trial=trial,
            arduino_time_s=t_s,
            raw=s,
            meta={"iti_s": iti_s}
        )

    # Cue start line => this is your best “trial onset” marker in the current log
    m = RE_CUE.search(s)
    if m:
        trial = int(m.group(1))
        label = m.group("label").strip()  # "Pair", "Stim only", "Tone only", etc.
        idx = m.group("idx")
        t_s = float(m.group("t"))
        meta = {"cue_label": label}
        if idx is not None:
            meta["cue_index"] = int(idx)

        # Treat cue start as TRIAL_ON for triggering purposes
        return Event(
            type="TRIAL_ON",
            wall_time_s=now,
            trial=trial,
            arduino_time_s=t_s,
            raw=s,
            meta=meta
        )

    # Outcome line
    m = RE_OUTCOME.search(s)
    if m:
        trial = int(m.group(1))
        label = m.group("label").strip()  # "small reward", "big reward", "small punish", "omission", etc.
        t_s = float(m.group("t"))
        return Event(
            type="OUTCOME",
            wall_time_s=now,
            trial=trial,
            arduino_time_s=t_s,
            raw=s,
            meta={"outcome_label": label}
        )

    # Timeout line
    m = RE_TIMEOUT.search(s)
    if m:
        trial = int(m.group(1))
        t_s = float(m.group("t"))
        return Event(
            type="TIMEOUT",
            wall_time_s=now,
            trial=trial,
            arduino_time_s=t_s,
            raw=s,
            meta={}
        )

    # Lick detected
    m = RE_LICK.search(s)
    if m:
        trial = int(m.group(1))
        t_s = float(m.group("t"))
        return Event(
            type="LICK",
            wall_time_s=now,
            trial=trial,
            arduino_time_s=t_s,
            raw=s,
            meta={}
        )

    # Simple ITI from OperantClamp_v1 (no trial number)
    m = RE_SIMPLE_ITI.search(s)
    if m:
        iti_s = float(m.group(1))
        return Event(
            type="ITI",
            wall_time_s=now,
            trial=None,
            arduino_time_s=None,
            raw=s,
            meta={"iti_s": iti_s},
        )

    # Lick from OperantClamp_v1 (no trial number)
    m = RE_SIMPLE_LICK.search(s)
    if m:
        t_s = float(m.group(1))
        return Event(
            type="LICK",
            wall_time_s=now,
            trial=None,
            arduino_time_s=t_s,
            raw=s,
            meta={},
        )

    # DA unclamp / clamp events from OperantClamp_v1
    m = RE_DA_UNCLAMP.search(s)
    if m:
        t_s = float(m.group("t"))
        count = int(m.group("count"))
        return Event(
            type="DA_UNCLAMP",
            wall_time_s=now,
            trial=None,
            arduino_time_s=t_s,
            raw=s,
            meta={"count": count},
        )

    m = RE_DA_CLAMP.search(s)
    if m:
        t_s = float(m.group("t"))
        count = int(m.group("count"))
        return Event(
            type="DA_CLAMP",
            wall_time_s=now,
            trial=None,
            arduino_time_s=t_s,
            raw=s,
            meta={"count": count},
        )

    m = RE_MANUAL_DA_UNCLAMP.search(s)
    if m:
        t_s = float(m.group("t"))
        count = int(m.group("count"))
        return Event(
            type="MANUAL_DA_UNCLAMP",
            wall_time_s=now,
            trial=None,
            arduino_time_s=t_s,
            raw=s,
            meta={"count": count},
        )

    m = RE_MANUAL_DA_CLAMP.search(s)
    if m:
        t_s = float(m.group("t"))
        count = int(m.group("count"))
        return Event(
            type="MANUAL_DA_CLAMP",
            wall_time_s=now,
            trial=None,
            arduino_time_s=t_s,
            raw=s,
            meta={"count": count},
        )

    # Task started (millis in ms, not divided by 1000)
    m = RE_TASK_START.search(s)
    if m:
        ms = int(m.group(1))
        return Event(
            type="TASK_START",
            wall_time_s=now,
            trial=None,
            arduino_time_s=ms / 1000.0,
            raw=s,
            meta={}
        )

    # Task ended (already in seconds)
    m = RE_TASK_END.search(s)
    if m:
        t_s = float(m.group(1))
        return Event(
            type="TASK_END",
            wall_time_s=now,
            trial=None,
            arduino_time_s=t_s,
            raw=s,
            meta={}
        )

    # Manual reward
    m = RE_MANUAL_REWARD.search(s)
    if m:
        timer_ms = int(m.group(1))
        size_ms = int(m.group(2))
        return Event(
            type="MANUAL_REWARD",
            wall_time_s=now,
            trial=None,
            arduino_time_s=timer_ms / 1000.0,
            raw=s,
            meta={"size_ms": size_ms}
        )

    # Manual punishment
    m = RE_MANUAL_PUNISH.search(s)
    if m:
        timer_ms = int(m.group(1))
        size_ms = int(m.group(2))
        return Event(
            type="MANUAL_PUNISH",
            wall_time_s=now,
            trial=None,
            arduino_time_s=timer_ms / 1000.0,
            raw=s,
            meta={"size_ms": size_ms}
        )

    # Not recognized
    return None


def cmd_pid_on(ramp_ms: Optional[int] = None) -> str:
    if ramp_ms is None:
        return "CMD PID_ON"
    return f"CMD PID_ON {int(ramp_ms)}"


def cmd_pid_off(ramp_ms: Optional[int] = None) -> str:
    if ramp_ms is None:
        return "CMD PID_OFF"
    return f"CMD PID_OFF {int(ramp_ms)}"


def cmd_set_target(value: float, dur: Optional[float] = None,
                   decay_mode: Optional[str] = None) -> str:
    """
    T0-relative df/F offset as a fraction (e.g., 0.50 means raw F = 1.50*T0_raw;
    0.0 = match T0). GUI applies new_target_dff = base + value*(100 + base);
    forwards to firmware as S<computed_dff>.
    If dur > 0 (seconds), the GUI holds that offset for dur then restores
    the previous absolute PID target.

    decay_mode in {"fopdt", "natural"} replaces the post-dur restore
    with an Arduino-side decay trajectory (Y command).  Channel for the decay
    model is selected by the GUI from sign(value): >0 -> excite_decay, <0 -> inhibit_decay.
    """
    if decay_mode is not None and decay_mode not in ("fopdt", "natural"):
        raise ValueError(f"decay_mode must be one of fopdt/natural, got {decay_mode!r}")
    base = f"CMD SET_TARGET {value:.4f}"
    if dur is not None and float(dur) > 0:
        base = f"{base} {float(dur):.6g}"
        if decay_mode is not None:
            base = f"{base} {decay_mode}"
    return base


def cmd_gui_set_target(dur_s: float) -> str:
    """
    Ask brainclamp_gui to compute the median DA value over the *next* dur_s seconds
    (based on its incoming DATA: stream) and use that as the PID target.
    """
    return f"CMD GUI_SET_TARGET {float(dur_s):.6g}"


def cmd_reset_baseline() -> str:
    return "CMD RESET_BASELINE"


def cmd_pulse_excite(amp: float, dur_s: float) -> str:
    """
    Excitatory pulse command.
      amp   : arbitrary amplitude units understood by brainclamp_gui/Arduino
      dur_s : duration in seconds
    """
    dur_ms = int(round(dur_s * 1000.0))
    return f"CMD PULSE_EXC {amp:.6g} {dur_ms}"


def cmd_pulse_inhibit(amp: float, dur_s: float) -> str:
    """
    Inhibitory pulse command.
      amp   : arbitrary amplitude units understood by brainclamp_gui/Arduino
      dur_s : duration in seconds
    """
    dur_ms = int(round(dur_s * 1000.0))
    return f"CMD PULSE_INH {amp:.6g} {dur_ms}"


_trial_response_window: Optional[int] = None


def apply_session_logic(ev: Event, udp: Optional[UdpSender]) -> None:
    """
    Customize this function per session/experiment.

    It is called for every parsed Event. Use the helpers above (cmd_pid_on,
    cmd_pid_off, cmd_set_target, cmd_pulse_excite, cmd_pulse_inhibit,
    cmd_reset_baseline) to decide which UDP commands to send to
    brainclamp_gui.

    Example patterns (uncomment and edit to taste):

        if ev.type == "TRIAL_ON":
            udp.send(cmd_pid_on())
            udp.send(cmd_set_target(0.5))   # +50% raw F over T0

        if ev.type == "OUTCOME" and ev.meta.get("outcome_label") == "big reward":
            udp.send(cmd_pulse_excite(amp=0.5, dur_s=0.5))

        if ev.type == "OUTCOME" and ev.meta.get("outcome_label") == "big punish":
            udp.send(cmd_pulse_inhibit(amp=0.5, dur_s=0.5))

        if ev.type == "TASK_END":
            udp.send(cmd_pid_off())

    """
    global _trial_response_window

    if udp is None:
        return

    # Turn on clamping when trial starts and open a response window for lick-triggered DA pulses.
    if ev.type == "TRIAL_ON":
        tone_delay_s = 0.4
        _trial_response_window = ev.trial
        send_protocol(udp, cmd_pid_on(ramp_ms=200))
        delay_send(udp, cmd_set_target(0.80, dur=0.25, decay_mode="natural"), tone_delay_s)   # +80% raw F over T0
        return

    # # Every lick inside the response window: drive a 100 ms phasic excitatory pulse.
    # if (
    #     ev.type == "LICK"
    #     and _trial_response_window is not None
    #     and ev.trial == _trial_response_window
    # ):
    #     send_protocol(udp, cmd_pulse_excite(amp=0.5, dur_s=0.1))
    #     return

    # Outcome (or timeout) closes the response window and shapes the post-outcome DA.
    if ev.type in ("OUTCOME", "TIMEOUT"):
        _trial_response_window = None

        if ev.meta.get("outcome_label") in ("big reward", "small reward"):
            send_protocol(udp, cmd_set_target(0.80))   # +80% raw F over T0
        elif ev.meta.get("outcome_label") == "omission":
            send_protocol(udp, cmd_set_target(-0.10, dur=0.25, decay_mode="natural"))

        delay_send(udp, cmd_pid_off(ramp_ms=500), 0.5)
        return
