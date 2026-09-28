"""Manual button layout, commands, and Arduino-reported toggle state."""
from dataclasses import dataclass
import re

CAPABILITY_LINE = "BEHAVIOR_CONTROLS 1"
COMMANDS = {
    "small_reward": "w", "large_reward": "W", "punishment": "p", "tone": "t",
    "blue_open": "b", "blue_close": "B", "blue_pattern": "f",
    "red_open": "r", "red_close": "R", "red_pattern": "F",
    "calibration": "c", "start_task": "s", "end_task": "x",
}
CONTROL_GROUPS = (
    ("Reward / sound", (("small_reward", "Small reward"),
                        ("large_reward", "Large reward"),
                        ("punishment", "Punishment (no tone)"), ("tone", "Tone"))),
    ("Lasers", (("blue_toggle", "Open blue"), ("red_toggle", "Open red"),
                ("blue_pattern", "Blue pattern"), ("red_pattern", "Red pattern"))),
    ("Session", (("calibration", "Water calibration"), ("task_toggle", "Start task"))),
)


@dataclass
class ManualState:
    task_running: bool = False
    blue_open: bool = False
    red_open: bool = False

    def consume(self, line):
        """Update only on board feedback, never optimistically on a GUI click."""
        snapshot = re.fullmatch(r"BEHAVIOR_STATE task=([01]) blue=([01]) red=([01])", line)
        if snapshot:
            self.task_running, self.blue_open, self.red_open = (v == "1" for v in snapshot.groups())
        elif line == "ACK START_TASK" or line.startswith("TASK STARTED AT"):
            self.task_running = True
            self.blue_open = self.red_open = False
        elif line.startswith("ACK END_TASK") or line.startswith("TASK ENDED AT"):
            self.task_running = self.blue_open = self.red_open = False
        elif line in ("ACK BLUE_OPEN", "ACK BLUE_PATTERN"):
            self.blue_open = True
        elif line in ("ACK RED_OPEN", "ACK RED_PATTERN"):
            self.red_open = True
        elif line in ("ACK BLUE_CLOSE", "DONE BLUE_PATTERN"):
            self.blue_open = False
        elif line in ("ACK RED_CLOSE", "DONE RED_PATTERN"):
            self.red_open = False
        else:
            return False
        return True

    def action_for(self, key):
        if key == "blue_toggle":
            return "blue_close" if self.blue_open else "blue_open"
        if key == "red_toggle":
            return "red_close" if self.red_open else "red_open"
        if key == "task_toggle":
            return "end_task" if self.task_running else "start_task"
        return key
