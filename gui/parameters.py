"""Edit explicitly marked Arduino settings without modifying source sketches.

This deliberately handles a small declaration/initializer syntax, not arbitrary
C++. Arduino CLI remains responsible for C++ type checking and compilation.
"""
from __future__ import annotations

import ast
from collections import Counter
from dataclasses import dataclass
import hashlib
from pathlib import Path
import re
import shutil


START = re.compile(r"^\s*//[^\n]*(?:user settings|GUI parameters begin)[^\n]*$", re.I | re.M)
END = re.compile(r"^\s*//[^\n]*(?:(?:params?|parameters?)\s+initial|GUI parameters end)[^\n]*$", re.I | re.M)
DECLARATION = re.compile(
    r"^[ \t]*(?P<type>(?:const\s+)?(?:unsigned\s+long(?:\s+int)?|unsigned\s+int|"
    r"unsigned\s+short|long(?:\s+int)?|short(?:\s+int)?|int|float|double|bool|boolean|byte|u?int(?:8|16|32)_t))"
    r"[ \t]+(?P<name>[A-Za-z_]\w*)[ \t]*(?:\[[ \t]*(?P<size>\d+)[ \t]*\])?"
    r"[ \t]*=[ \t]*(?P<value>[^;\r\n]+?)[ \t]*;[ \t\r\n]*$"
)
COMMENTS_AND_STRINGS = re.compile(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'')
NUMBER = re.compile(r"(?:0[xX][0-9a-fA-F]+|(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)(?:[uUlL]+|[fF])?")


@dataclass(frozen=True)
class Parameter:
    name: str
    ctype: str
    value: str
    start: int
    end: int
    group: str = "Behavior settings"
    comment: str = ""
    size: int | None = None

    @property
    def boolean(self):
        return self.ctype.removeprefix("const ") in ("bool", "boolean")


def source_hash(data):
    return hashlib.sha256(data).hexdigest()


def parameters_in(source):
    """Return unambiguous single-line globals inside the marked settings block."""
    begin = START.search(source)
    end = END.search(source, begin.end()) if begin else None
    if not begin or not end:
        return []
    # Preserve offsets while hiding comments and strings from the recognizer.
    masked = COMMENTS_AND_STRINGS.sub(lambda m: re.sub(r"[^\r\n]", " ", m[0]), source)
    offset = begin.end()
    depth = conditional = 0
    group = "Behavior settings"
    found = []
    for line in masked[offset:end.start()].splitlines(keepends=True):
        original = source[offset:offset + len(line)]
        stripped = line.strip()
        if re.match(r"#\s*(if|ifdef|ifndef)\b", stripped):
            conditional += 1
        elif re.match(r"#\s*endif\b", stripped):
            conditional = max(0, conditional - 1)
        if not conditional and depth == 0:
            heading = re.match(r"\s*//\s*(.*(?:params?|parameters?|settings).*)", original, re.I)
            if heading:
                group = heading[1].strip(" */\r\n")
            match = DECLARATION.fullmatch(line)
            if match:
                start, finish = offset + match.start("value"), offset + match.end("value")
                value = source[start:finish]
                # A string/comment inside the initializer is outside supported syntax.
                if value == match["value"]:
                    comment = original.split("//", 1)[1].strip() if "//" in original else ""
                    parameter = Parameter(match["name"], " ".join(match["type"].split()), value,
                                          start, finish, group, comment,
                                          int(match["size"]) if match["size"] else None)
                    try:
                        validate_value(parameter, value)
                    except ValueError:
                        # In particular, never treat "int a = 1, b = 2;" as
                        # one initializer: replacing it would delete b.
                        pass
                    else:
                        found.append(parameter)
        if not conditional:
            depth = max(0, depth + line.count("{") - line.count("}"))
        offset += len(line)
    counts = Counter(p.name for p in found)
    return [p for p in found if counts[p.name] == 1]


def _expression(value):
    """Accept numbers, names, parentheses and basic arithmetic; never eval input."""
    if not value or len(value) > 200:
        raise ValueError("Enter a value or a short arithmetic expression.")
    tokens = []
    pos = 0
    while pos < len(value):
        if value[pos].isspace():
            pos += 1
            continue
        number = NUMBER.match(value, pos)
        name = re.match(r"[A-Za-z_]\w*", value[pos:])
        if number:
            token = number[0]
            # Strip C++ suffixes for syntax checking; keep original text for Arduino.
            normalized = re.sub(r"[uUlL]+$", "", token)
            if not normalized.lower().startswith("0x"):
                normalized = re.sub(r"[fF]$", "", normalized)
            tokens.append(normalized)
            pos = number.end()
        elif name:
            tokens.append(name[0])
            pos += len(name[0])
        elif value[pos] in "+-*/%()":
            tokens.append(value[pos])
            pos += 1
        else:
            raise ValueError("Use numbers, parameter names, + - * / %, and parentheses only.")
    if "//" in value or "/*" in value or "**" in value:
        raise ValueError("Comments and exponent operators are not supported.")
    try:
        tree = ast.parse(" ".join(tokens), mode="eval")
    except SyntaxError as exc:
        raise ValueError("Invalid numeric expression.") from exc
    allowed = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Constant, ast.Name, ast.Load,
               ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Mod, ast.UAdd, ast.USub)
    if any(not isinstance(node, allowed) for node in ast.walk(tree)):
        raise ValueError("Only numeric arithmetic expressions are supported.")


def validate_value(parameter, value):
    value = value.strip()
    try:
        if parameter.size is not None:
            if not (value.startswith("{") and value.endswith("}")):
                raise ValueError("Enter an array in braces, for example {1, 100}.")
            items = [item.strip() for item in value[1:-1].split(",")]
            if len(items) != parameter.size:
                raise ValueError(f"Enter exactly {parameter.size} array values.")
        else:
            items = [value]
        for item in items:
            if parameter.boolean:
                if item not in ("true", "false"):
                    raise ValueError("Choose true or false.")
            else:
                _expression(item)
                # Reject obvious negative unsigned literals; expression types are
                # checked by the compiler, since int widths vary between boards.
                if ("unsigned" in parameter.ctype or "uint" in parameter.ctype or parameter.ctype == "byte") and re.fullmatch(r"-\s*\d+[uUlL]*", item):
                    raise ValueError("Unsigned values cannot be negative.")
        if parameter.name == "UnitRewardSize":
            if not re.fullmatch(r"[0-9]+", value) or not 1 <= int(value) <= 0xFFFFFFFF:
                raise ValueError("Enter a positive whole number of milliseconds (1–4294967295).")
    except ValueError as exc:
        raise ValueError(f"{parameter.name}: {exc}") from exc
    return value


def apply_overrides(data, overrides, expected_hash):
    """Patch initializer spans only; reject stale snapshots and unknown names."""
    if not isinstance(overrides, dict) or any(not isinstance(k, str) or not isinstance(v, str)
                                             for k, v in overrides.items()):
        raise ValueError("Parameter overrides must map parameter names to text values.")
    if not overrides:
        return data
    if source_hash(data) != expected_hash:
        raise ValueError("The Arduino source has changed since these parameters were loaded. "
                         "Reload from sketch, review the new values, and reapply your edits.")
    source = data.decode("utf-8")
    parameters = {p.name: p for p in parameters_in(source)}
    unknown = set(overrides) - set(parameters)
    if unknown:
        raise ValueError("Parameters not editable in this sketch: " + ", ".join(sorted(unknown)))
    replacements = [(parameters[name], validate_value(parameters[name], value))
                    for name, value in overrides.items()]
    for parameter, value in sorted(replacements, key=lambda item: item[0].start, reverse=True):
        source = source[:parameter.start] + value + source[parameter.end:]
    return source.encode("utf-8")


def prepare_sketch(sketch, destination, overrides, expected_hash):
    """Copy the whole sketch, then modify only the main .ino in that copy."""
    target = Path(destination) / sketch.name
    shutil.copytree(sketch, target)
    main = target / (sketch.name + ".ino")
    data = main.read_bytes()
    updated = apply_overrides(data, overrides, expected_hash)
    if updated != data:
        main.write_bytes(updated)
    return target


class ParameterEditor:
    """Pinned calibration field plus a scrollable editor for the remaining values."""
    def __init__(self, parent, reload_command):
        import tkinter as tk
        from tkinter import ttk
        self.tk, self.ttk = tk, ttk
        self.frame = ttk.LabelFrame(parent, text="Arduino parameters", padding=10)
        self.frame.columnconfigure(0, weight=1)
        self.frame.rowconfigure(3, weight=1)
        self.parameters = []
        self.values = {}
        self.widgets = []
        self.hash = ""
        self.enabled = True
        self.source = b""
        calibration = ttk.Frame(self.frame)
        calibration.grid(row=0, column=0, sticky="ew")
        calibration.columnconfigure(1, weight=1)
        ttk.Label(calibration, text="UnitRewardSize", font=("Helvetica", 11, "bold")).grid(row=0, column=0, sticky="w", padx=(0, 8))
        self.unit = tk.StringVar()
        self.unit_entry = ttk.Entry(calibration, textvariable=self.unit, width=10)
        self.unit_entry.grid(row=0, column=1, sticky="ew")
        ttk.Label(calibration, text="ms").grid(row=0, column=2, padx=(5, 0))
        ttk.Label(self.frame, text="Daily valve calibration · applies on upload", style="Hint.TLabel").grid(row=1, column=0, sticky="w", pady=(4, 6))
        self.reload_button = ttk.Button(self.frame, text="Reload from sketch / reset edits", command=reload_command)
        self.reload_button.grid(row=2, column=0, sticky="ew", pady=(0, 8))
        scroller = ttk.Frame(self.frame)
        scroller.grid(row=3, column=0, sticky="nsew")
        scroller.rowconfigure(0, weight=1)
        scroller.columnconfigure(0, weight=1)
        self.canvas = tk.Canvas(scroller, height=255, width=340, highlightthickness=0)
        self.canvas.grid(row=0, column=0, sticky="nsew")
        scrollbar = ttk.Scrollbar(scroller, orient="vertical", command=self.canvas.yview)
        scrollbar.grid(row=0, column=1, sticky="ns")
        self.canvas.configure(yscrollcommand=scrollbar.set)
        self.body = ttk.Frame(self.canvas)
        self.body.columnconfigure(1, weight=1)
        self.window = self.canvas.create_window((0, 0), window=self.body, anchor="nw")
        self.body.bind("<Configure>", lambda event: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.bind("<Configure>", lambda event: self.canvas.itemconfigure(self.window, width=event.width))
        self.note = tk.StringVar()
        ttk.Label(self.frame, textvariable=self.note, style="Hint.TLabel", wraplength=350).grid(row=4, column=0, sticky="w", pady=(6, 0))
        self.clear("Select an Arduino sketch to load its parameters.")

    def clear(self, message):
        for widget in self.body.winfo_children():
            widget.destroy()
        self.parameters, self.values, self.widgets = [], {}, []
        self.hash, self.source = "", b""
        self.unit.set("Not defined")
        self.unit_entry.configure(state="disabled")
        self.note.set(message)

    def load(self, main, overrides=None, expected_hash=""):
        data = Path(main).read_bytes()
        overrides = overrides or {}
        apply_overrides(data, overrides, expected_hash)  # validate profile before changing the editor
        parameters = parameters_in(data.decode("utf-8"))
        self.clear("")
        self.source, self.hash, self.parameters = data, source_hash(data), parameters
        row = 0
        previous_group = None
        for parameter in parameters:
            value = overrides.get(parameter.name, parameter.value)
            if parameter.name == "UnitRewardSize":
                self.unit.set(value)
                self.values[parameter.name] = self.unit
                continue
            if parameter.group != previous_group:
                self.ttk.Label(self.body, text=parameter.group, font=("Helvetica", 10, "bold"),
                               wraplength=325).grid(row=row, column=0, columnspan=2, sticky="w", pady=(10, 5))
                previous_group = parameter.group
                row += 1
            variable = self.tk.StringVar(value=value)
            self.values[parameter.name] = variable
            label = parameter.name + (f" [{parameter.size}]" if parameter.size is not None else "")
            self.ttk.Label(self.body, text=label, wraplength=170).grid(row=row, column=0, sticky="w", padx=(0, 6), pady=3)
            boolean = parameter.boolean and parameter.size is None
            widget = (self.ttk.Combobox(self.body, values=("true", "false"), state="readonly", textvariable=variable, width=17)
                      if boolean else self.ttk.Entry(self.body, textvariable=variable, width=17))
            widget.grid(row=row, column=1, sticky="ew", pady=3)
            self.widgets.append((widget, boolean))
            row += 1
            if parameter.comment:
                self.ttk.Label(self.body, text=parameter.comment, style="Hint.TLabel", wraplength=325).grid(
                    row=row, column=0, columnspan=2, sticky="w", pady=(0, 5))
                row += 1
        self.note.set(f"{len(parameters)} parameters · expressions are preserved."
                      if parameters else "No marked user settings found. See the setup guide to mark editable parameters.")
        self.canvas.yview_moveto(0)
        self.set_enabled(self.enabled)
        self._bind_scrolling(self.canvas)

    def _bind_scrolling(self, widget):
        widget.bind("<MouseWheel>", self._scroll)
        widget.bind("<Button-4>", self._scroll)
        widget.bind("<Button-5>", self._scroll)
        for child in widget.winfo_children():
            self._bind_scrolling(child)

    def _scroll(self, event):
        if getattr(event, "num", None) in (4, 5):
            amount = -1 if event.num == 4 else 1
        else:
            delta = event.delta
            amount = -int(delta / 120) if abs(delta) >= 120 else (-1 if delta > 0 else 1)
        self.canvas.yview_scroll(amount, "units")
        return "break"

    def set_enabled(self, enabled):
        self.enabled = enabled
        self.reload_button.configure(state="normal" if enabled else "disabled")
        self.unit_entry.configure(state="normal" if enabled and "UnitRewardSize" in self.values else "disabled")
        for widget, boolean in self.widgets:
            widget.configure(state=("readonly" if boolean else "normal") if enabled else "disabled")

    def overrides(self):
        result = {}
        for parameter in self.parameters:
            value = self.values[parameter.name].get().strip()
            if value != parameter.value:
                result[parameter.name] = validate_value(parameter, value)
        return result
