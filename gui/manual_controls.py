"""Version 1 manual commands shared by the two DA-clamp sketches."""

CAPABILITY_LINE = "BEHAVIOR_CONTROLS 1"
# Dedicated single-byte commands avoid the incompatible legacy digit mappings.
CONTROL_GROUPS = (
    ("Reward / sound", (("small_reward", "Small reward", "w"),
                        ("large_reward", "Large reward", "W"),
                        ("punishment", "Punishment (no tone)", "p"),
                        ("tone", "Tone", "t"))),
    ("Blue laser", (("blue_open", "Open blue", "b"),
                    ("blue_close", "Close blue", "B"),
                    ("blue_pattern", "Blue pattern", "f"))),
    ("Red laser", (("red_open", "Open red", "r"),
                   ("red_close", "Close red", "R"),
                   ("red_pattern", "Red pattern", "F"))),
    ("Session", (("calibration", "Water calibration", "c"),
                 ("start_task", "Start task", "s"),
                 ("end_task", "End task (keep counts)", "x"))),
)
COMMANDS = {key: command for _, group in CONTROL_GROUPS for key, _, command in group}
