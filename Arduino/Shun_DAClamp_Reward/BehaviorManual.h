#ifndef NEURODAP_BEHAVIOR_MANUAL_H
#define NEURODAP_BEHAVIOR_MANUAL_H
// Keep all standalone copies of this header identical.
// Each sketch is self-contained for Arduino IDE and temporary GUI builds.
#ifndef BEHAVIOR_LEGACY_SECOND_TONE
#define BEHAVIOR_LEGACY_SECOND_TONE (!BEHAVIOR_LEGACY_REWARD_BIG)
#endif
#ifndef BEHAVIOR_SHUTTER_OPEN
#define BEHAVIOR_SHUTTER_OPEN LOW
#define BEHAVIOR_SHUTTER_CLOSED HIGH
#endif
#ifndef BEHAVIOR_TASK_VALID
#define BEHAVIOR_TASK_VALID true
#endif
#ifndef BEHAVIOR_HAS_WATER
#define BEHAVIOR_HAS_WATER 1
#endif
#ifndef BEHAVIOR_HAS_AIR
#define BEHAVIOR_HAS_AIR 1
#endif
#ifndef BEHAVIOR_HAS_TONE
#define BEHAVIOR_HAS_TONE 1
#endif
#ifndef BEHAVIOR_HAS_TONE_COPIES
#define BEHAVIOR_HAS_TONE_COPIES 1
#endif
#ifndef BEHAVIOR_HAS_BLUE
#define BEHAVIOR_HAS_BLUE 1
#endif
#ifndef BEHAVIOR_HAS_RED
#define BEHAVIOR_HAS_RED 1
#endif
#ifndef BEHAVIOR_HAS_TASK
#define BEHAVIOR_HAS_TASK 1
#endif
void behaviorStartTask();
void behaviorStopTask();

struct BehaviorPulse {
  bool active;
  unsigned long started;
  unsigned long duration;
};
struct BehaviorPattern {
  bool active;
  bool on;
  unsigned long changed;
  unsigned int remaining;
};
BehaviorPulse behaviorWater = {false, 0, 0};
BehaviorPulse behaviorAir = {false, 0, 0};
BehaviorPattern behaviorBlue = {false, false, 0, 0};
BehaviorPattern behaviorRed = {false, false, 0, 0};
bool behaviorToneActive = false;
unsigned long behaviorToneStarted = 0;
bool behaviorCalibrating = false;
unsigned int behaviorCalibrationDelivered = 0;
unsigned long behaviorCalibrationChanged = 0;
bool behaviorCalibrationWaterOn = false;
#if BEHAVIOR_HAS_WATER
byte behaviorCalibrationValve = WaterSpout2;
byte behaviorCalibrationCopy = WaterSpout2_copy;
unsigned long behaviorCalibrationSize = 0;
unsigned int behaviorCalibrationRepeats = 0;
unsigned long behaviorCalibrationInterval = 0;
#endif

void behaviorReportState() {
  Serial.print("BEHAVIOR_STATE task="); Serial.print(state != Idle ? 1 : 0);
  Serial.print(" blue=");
#if BEHAVIOR_HAS_BLUE
  Serial.print((digitalRead(ShutterBlue) == BEHAVIOR_SHUTTER_OPEN || behaviorBlue.active) ? 1 : 0);
#else
  Serial.print(0);
#endif
  Serial.print(" red=");
#if BEHAVIOR_HAS_RED
  Serial.print((digitalRead(ShutterRed) == BEHAVIOR_SHUTTER_OPEN || behaviorRed.active) ? 1 : 0);
#else
  Serial.print(0);
#endif
  Serial.print(" calibration="); Serial.println(behaviorCalibrating ? 1 : 0);
}

void behaviorCapabilities() {
  Serial.println("BEHAVIOR_CONTROLS 1");
}

void behaviorStopOutputs() {
  // Cancel actuator work, not recorded session or calibration counts.
  behaviorWater.active = false;
  behaviorAir.active = false;
  behaviorBlue.active = false;
  behaviorRed.active = false;
  behaviorToneActive = false;
  behaviorCalibrating = false;
  behaviorCalibrationWaterOn = false;
#if BEHAVIOR_HAS_WATER
  LeftOutcomeButton = 0;
  RightOutcomeButton = 0;
  digitalWrite(WaterSpout, LOW);
  digitalWrite(WaterSpout_copy, LOW);
  digitalWrite(WaterSpout2, LOW);
  digitalWrite(WaterSpout2_copy, LOW);
#endif
#if BEHAVIOR_HAS_AIR
  digitalWrite(Airpuff, LOW);
  digitalWrite(Airpuff_copy, LOW);
#endif
#if BEHAVIOR_HAS_BLUE
  digitalWrite(ShutterBlue, BEHAVIOR_SHUTTER_CLOSED);
#endif
#if BEHAVIOR_HAS_RED
  digitalWrite(ShutterRed, BEHAVIOR_SHUTTER_CLOSED);
#endif
#if BEHAVIOR_HAS_TONE
  noTone(Speaker);
#if BEHAVIOR_HAS_TONE_COPIES
  digitalWrite(SpeakerLeft_copy, LOW);
  digitalWrite(SpeakerRight_copy, LOW);
#endif
#endif
}

void behaviorBeginTone() {
#if BEHAVIOR_HAS_TONE
  behaviorToneActive = true;
  behaviorToneStarted = millis();
  tone(Speaker, LeftCueFreq);
#if BEHAVIOR_HAS_TONE_COPIES
  digitalWrite(SpeakerLeft_copy, HIGH);
#endif
#endif
}

bool behaviorStartPattern(BehaviorPattern &pattern, byte pin, unsigned long pulse,
                          unsigned long period, unsigned int count) {
  if (pulse == 0 || period < pulse || count == 0) {
    Serial.println("ERR PATTERN_INVALID check pulse width, period, and count");
    return false;
  }
  pattern.active = true;
  pattern.on = true;
  pattern.changed = millis();
  pattern.remaining = count;
  digitalWrite(pin, BEHAVIOR_SHUTTER_OPEN);
  return true;
}

void behaviorUpdatePattern(BehaviorPattern &pattern, byte pin,
                           unsigned long pulse, unsigned long period) {
  if (!pattern.active) return;
  unsigned long duration = pattern.on ? pulse : period - pulse;
  if (millis() - pattern.changed < duration) return;
  pattern.changed = millis();
  if (pattern.on) {
    digitalWrite(pin, BEHAVIOR_SHUTTER_CLOSED);
    pattern.on = false;
    if (--pattern.remaining == 0) {
      pattern.active = false;
      Serial.println(&pattern == &behaviorBlue ? "DONE BLUE_PATTERN" : "DONE RED_PATTERN");
    }
  } else {
    digitalWrite(pin, BEHAVIOR_SHUTTER_OPEN);
    pattern.on = true;
  }
}

void behaviorUpdateManual() {
#if BEHAVIOR_HAS_WATER
  if (behaviorWater.active && millis() - behaviorWater.started >= behaviorWater.duration) {
    digitalWrite(WaterSpout2, LOW);
    digitalWrite(WaterSpout2_copy, LOW);
    behaviorWater.active = false;
    Serial.print("Manual reward\t");
    Serial.print(behaviorWater.started);
    Serial.print("\t");
    Serial.println(behaviorWater.duration);
  }
#endif
#if BEHAVIOR_HAS_AIR
  if (behaviorAir.active && millis() - behaviorAir.started >= behaviorAir.duration) {
    digitalWrite(Airpuff, LOW);
    digitalWrite(Airpuff_copy, LOW);
    behaviorAir.active = false;
    Serial.print("Manual punishment\t");
    Serial.print(behaviorAir.started);
    Serial.print("\t");
    Serial.println(behaviorAir.duration);
  }
#endif
#if BEHAVIOR_HAS_TONE
  if (behaviorToneActive && millis() - behaviorToneStarted >= ShortToneDuration) {
    noTone(Speaker);
#if BEHAVIOR_HAS_TONE_COPIES
    digitalWrite(SpeakerLeft_copy, LOW);
#endif
    behaviorToneActive = false;
  }
#endif
#if BEHAVIOR_HAS_BLUE
  behaviorUpdatePattern(behaviorBlue, ShutterBlue, BluePatternPulseMs, BluePatternPeriodMs);
#endif
#if BEHAVIOR_HAS_RED
  behaviorUpdatePattern(behaviorRed, ShutterRed, RedPatternPulseMs, RedPatternPeriodMs);
#endif
#if BEHAVIOR_HAS_WATER
  if (behaviorCalibrating) {
    unsigned long duration = behaviorCalibrationWaterOn ? behaviorCalibrationSize : behaviorCalibrationInterval;
    if (millis() - behaviorCalibrationChanged >= duration) {
      behaviorCalibrationChanged = millis();
      if (behaviorCalibrationWaterOn) {
        digitalWrite(behaviorCalibrationValve, LOW);
        digitalWrite(behaviorCalibrationCopy, LOW);
        behaviorCalibrationWaterOn = false;
        ++behaviorCalibrationDelivered;
        if (behaviorCalibrationDelivered >= behaviorCalibrationRepeats) {
          behaviorCalibrating = false;
          Serial.print("DONE CALIBRATION deliveries=");
          Serial.println(behaviorCalibrationDelivered);
        }
      } else {
        digitalWrite(behaviorCalibrationValve, HIGH);
        digitalWrite(behaviorCalibrationCopy, HIGH);
        behaviorCalibrationWaterOn = true;
      }
    }
  }
#endif
}

void behaviorHandleCommand(char command) {
  if (command == '\r' || command == '\n' || command == ' ') return;
  if (command == '?') { behaviorCapabilities(); behaviorReportState(); return; }
#ifdef BEHAVIOR_CUSTOM_LEGACY
  if (behaviorHandleLegacy(command)) return;
#endif
  // Preserve legacy digit meanings; the GUI uses the unambiguous letters.
  if (command == '1') command = BEHAVIOR_LEGACY_REWARD_BIG ? 'W' : 'w';
  if (command == '2') command = BEHAVIOR_LEGACY_SECOND_TONE ? 't' : 'p';
  if (command == '3') command = 'b';
  if (command == '4') command = 'r';
  if (command == '5') command = 'f';
  if (command == '6') command = 'F';
  if (command == '7') command = 'c';
  if (command == '8') command = 's';
  if (command == '9') command = 'x';
  if ((!BEHAVIOR_HAS_WATER && (command == 'w' || command == 'W' || command == 'c')) ||
      (!BEHAVIOR_HAS_AIR && command == 'p') || (!BEHAVIOR_HAS_TONE && command == 't') ||
      (!BEHAVIOR_HAS_BLUE && (command == 'b' || command == 'B' || command == 'f')) ||
      (!BEHAVIOR_HAS_RED && (command == 'r' || command == 'R' || command == 'F')) ||
      (!BEHAVIOR_HAS_TASK && command == 's')) {
    Serial.println("ERR UNSUPPORTED this sketch does not define that output/task");
    return;
  }
  if (command == 'x') {
    behaviorStopTask();
    Serial.println("ACK END_TASK counts_preserved");
    return;
  }
  if (command == 'B' || command == 'R') {
#if BEHAVIOR_HAS_BLUE
    if (command == 'B') { behaviorBlue.active = false; digitalWrite(ShutterBlue, BEHAVIOR_SHUTTER_CLOSED); }
#endif
#if BEHAVIOR_HAS_RED
    if (command == 'R') { behaviorRed.active = false; digitalWrite(ShutterRed, BEHAVIOR_SHUTTER_CLOSED); }
#endif
    Serial.println(command == 'B' ? "ACK BLUE_CLOSE" : "ACK RED_CLOSE");
    return;
  }
  if (command != 'w' && command != 'W' && command != 'p' && command != 't' &&
      command != 'b' && command != 'r' && command != 'f' && command != 'F' &&
      command != 'c' && command != 's') {
    Serial.println("ERR UNKNOWN_COMMAND");
    return;
  }
  // Prevent manual stimuli and calibration from competing with trial outputs.
  if (state != Idle || behaviorCalibrating) {
    Serial.println("ERR BUSY stop the task/calibration before manual stimuli or restart");
    return;
  }
  if (command == 's') {
    if (!BEHAVIOR_TASK_VALID) {
      Serial.println("ERR TASK_CONFIG blue task stimulation conflicts with paAIP2");
      return;
    }
    behaviorStopOutputs();
    behaviorStartTask();
    Serial.println("ACK START_TASK");
    return;
  }
  if (command == 'w' || command == 'W') {
#if BEHAVIOR_HAS_WATER
    if (behaviorWater.active) { Serial.println("ERR WATER_BUSY"); return; }
    behaviorWater.active = true;
    behaviorWater.started = millis();
    behaviorWater.duration = command == 'w' ? SmallRewardSize : BigRewardSize;
    ++ManualPositiveNum;
    digitalWrite(WaterSpout2, HIGH);
    digitalWrite(WaterSpout2_copy, HIGH);
    Serial.println(command == 'w' ? "ACK SMALL_REWARD" : "ACK LARGE_REWARD");
#endif
  } else if (command == 'p') {
#if BEHAVIOR_HAS_AIR
    if (behaviorAir.active) { Serial.println("ERR AIR_BUSY"); return; }
    behaviorAir.active = true;
    behaviorAir.started = millis();
    behaviorAir.duration = SmallPunishSize;
    ++ManualNegativeNum;
    digitalWrite(Airpuff, HIGH);
    digitalWrite(Airpuff_copy, HIGH);
    Serial.println("ACK PUNISHMENT");
#endif
  } else if (command == 't') {
    behaviorBeginTone();
    Serial.println("ACK TONE");
  } else if (command == 'b' || command == 'r') {
#if BEHAVIOR_HAS_BLUE
    if (command == 'b') { behaviorBlue.active = false; digitalWrite(ShutterBlue, BEHAVIOR_SHUTTER_OPEN); }
#endif
#if BEHAVIOR_HAS_RED
    if (command == 'r') { behaviorRed.active = false; digitalWrite(ShutterRed, BEHAVIOR_SHUTTER_OPEN); }
#endif
    Serial.println(command == 'b' ? "ACK BLUE_OPEN" : "ACK RED_OPEN");
  } else if (command == 'f') {
#if BEHAVIOR_HAS_BLUE
    if (behaviorStartPattern(behaviorBlue, ShutterBlue, BluePatternPulseMs, BluePatternPeriodMs, BluePatternPulseCount))
      Serial.println("ACK BLUE_PATTERN");
#endif
  } else if (command == 'F') {
#if BEHAVIOR_HAS_RED
    if (behaviorStartPattern(behaviorRed, ShutterRed, RedPatternPulseMs, RedPatternPeriodMs, RedPatternPulseCount))
      Serial.println("ACK RED_PATTERN");
#endif
  } else if (command == 'c') {
#if BEHAVIOR_HAS_WATER
    if (behaviorWater.active || CalibrationRepeats == 0 || UnitRewardSize == 0) {
      Serial.println("ERR CALIBRATION_BUSY_OR_INVALID");
      return;
    }
    behaviorCalibrationValve = WaterSpout2;
    behaviorCalibrationCopy = WaterSpout2_copy;
    behaviorCalibrationSize = UnitRewardSize;
    behaviorCalibrationRepeats = CalibrationRepeats;
    behaviorCalibrationInterval = CalibrationIntervalMs;
    behaviorCalibrating = true;
    behaviorCalibrationDelivered = 0; // reset only when explicitly starting a new calibration
    behaviorCalibrationWaterOn = false;
    behaviorCalibrationChanged = millis();
    Serial.println("ACK CALIBRATION");
#endif
  }
}

void behaviorReadSerial() {
  // Bound work per loop so output timers still run under continuous input.
  for (byte n = 0; n < 16 && Serial.available() > 0; ++n) {
    behaviorHandleCommand((char)Serial.read());
  }
}
// Used by legacy task tones, which otherwise block serial stop handling.
bool behaviorWait(unsigned long duration) {
  unsigned long began = millis();
  while (true) {
    if (Serial.available() > 0) {
      char command = (char)Serial.read();
      if (command == 'x' || command == '9' || command == '?') behaviorHandleCommand(command);
      else if (command != '\r' && command != '\n' && command != ' ')
        Serial.println("ERR BUSY task tone in progress");
    }
    if (state == Idle) return false;
    if (millis() - began >= duration) return true;
    delay(1);
  }
}
#endif
