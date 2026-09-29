#ifndef NEURODAP_BEHAVIOR_MANUAL_H
#define NEURODAP_BEHAVIOR_MANUAL_H
// Keep the copies in Shun_DAClamp_Reward and Shun_DAClamp_Random identical.
// Each sketch is self-contained for Arduino IDE and temporary GUI builds.
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

void behaviorReportState() {
  Serial.print("BEHAVIOR_STATE task="); Serial.print(state != Idle ? 1 : 0);
  Serial.print(" blue="); Serial.print((digitalRead(ShutterBlue) == LOW || behaviorBlue.active) ? 1 : 0);
  Serial.print(" red="); Serial.print((digitalRead(ShutterRed) == LOW || behaviorRed.active) ? 1 : 0);
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
  LeftOutcomeButton = 0;
  RightOutcomeButton = 0;
  digitalWrite(WaterSpout, LOW);
  digitalWrite(WaterSpout_copy, LOW);
  digitalWrite(WaterSpout2, LOW);
  digitalWrite(WaterSpout2_copy, LOW);
  digitalWrite(Airpuff, LOW);
  digitalWrite(Airpuff_copy, LOW);
  digitalWrite(ShutterBlue, HIGH);
  digitalWrite(ShutterRed, HIGH);
  noTone(Speaker);
  digitalWrite(SpeakerLeft_copy, LOW);
  digitalWrite(SpeakerRight_copy, LOW);
}

void behaviorBeginTone() {
  behaviorToneActive = true;
  behaviorToneStarted = millis();
  tone(Speaker, LeftCueFreq);
  digitalWrite(SpeakerLeft_copy, HIGH);
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
  digitalWrite(pin, LOW);
  return true;
}

void behaviorUpdatePattern(BehaviorPattern &pattern, byte pin,
                           unsigned long pulse, unsigned long period) {
  if (!pattern.active) return;
  unsigned long duration = pattern.on ? pulse : period - pulse;
  if (millis() - pattern.changed < duration) return;
  pattern.changed = millis();
  if (pattern.on) {
    digitalWrite(pin, HIGH);
    pattern.on = false;
    if (--pattern.remaining == 0) {
      pattern.active = false;
      Serial.println(pin == ShutterBlue ? "DONE BLUE_PATTERN" : "DONE RED_PATTERN");
    }
  } else {
    digitalWrite(pin, LOW);
    pattern.on = true;
  }
}

void behaviorUpdateManual() {
  if (behaviorWater.active && millis() - behaviorWater.started >= behaviorWater.duration) {
    digitalWrite(WaterSpout2, LOW);
    digitalWrite(WaterSpout2_copy, LOW);
    behaviorWater.active = false;
    Serial.print("Manual reward\t");
    Serial.print(behaviorWater.started);
    Serial.print("\t");
    Serial.println(behaviorWater.duration);
  }
  if (behaviorAir.active && millis() - behaviorAir.started >= behaviorAir.duration) {
    digitalWrite(Airpuff, LOW);
    digitalWrite(Airpuff_copy, LOW);
    behaviorAir.active = false;
    Serial.print("Manual punishment\t");
    Serial.print(behaviorAir.started);
    Serial.print("\t");
    Serial.println(behaviorAir.duration);
  }
  if (behaviorToneActive && millis() - behaviorToneStarted >= ShortToneDuration) {
    noTone(Speaker);
    digitalWrite(SpeakerLeft_copy, LOW);
    behaviorToneActive = false;
  }
  behaviorUpdatePattern(behaviorBlue, ShutterBlue, BluePatternPulseMs, BluePatternPeriodMs);
  behaviorUpdatePattern(behaviorRed, ShutterRed, RedPatternPulseMs, RedPatternPeriodMs);
  if (behaviorCalibrating) {
    unsigned long duration = behaviorCalibrationWaterOn ? UnitRewardSize : CalibrationIntervalMs;
    if (millis() - behaviorCalibrationChanged >= duration) {
      behaviorCalibrationChanged = millis();
      if (behaviorCalibrationWaterOn) {
        digitalWrite(WaterSpout2, LOW);
        digitalWrite(WaterSpout2_copy, LOW);
        behaviorCalibrationWaterOn = false;
        ++behaviorCalibrationDelivered;
        if (behaviorCalibrationDelivered >= CalibrationRepeats) {
          behaviorCalibrating = false;
          Serial.print("DONE CALIBRATION deliveries=");
          Serial.println(behaviorCalibrationDelivered);
        }
      } else {
        digitalWrite(WaterSpout2, HIGH);
        digitalWrite(WaterSpout2_copy, HIGH);
        behaviorCalibrationWaterOn = true;
      }
    }
  }
}

void behaviorHandleCommand(char command) {
  if (command == '\r' || command == '\n' || command == ' ') return;
  if (command == '?') { behaviorCapabilities(); behaviorReportState(); return; }
  // Preserve legacy digit meanings; the GUI uses the unambiguous letters.
  if (command == '1') command = BEHAVIOR_LEGACY_REWARD_BIG ? 'W' : 'w';
  if (command == '2') command = BEHAVIOR_LEGACY_REWARD_BIG ? 'p' : 't';
  if (command == '3') command = 'b';
  if (command == '4') command = 'r';
  if (command == '5') command = 'f';
  if (command == '6') command = 'F';
  if (command == '7') command = 'c';
  if (command == '8') command = 's';
  if (command == '9') command = 'x';
  if (command == 'x') {
    behaviorStopTask();
    Serial.println("ACK END_TASK counts_preserved");
    return;
  }
  if (command == 'B' || command == 'R') {
    if (command == 'B') { behaviorBlue.active = false; digitalWrite(ShutterBlue, HIGH); }
    else { behaviorRed.active = false; digitalWrite(ShutterRed, HIGH); }
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
    behaviorStopOutputs();
    behaviorStartTask();
    Serial.println("ACK START_TASK");
    return;
  }
  if (command == 'w' || command == 'W') {
    if (behaviorWater.active) { Serial.println("ERR WATER_BUSY"); return; }
    behaviorWater.active = true;
    behaviorWater.started = millis();
    behaviorWater.duration = command == 'w' ? SmallRewardSize : BigRewardSize;
    ++ManualPositiveNum;
    digitalWrite(WaterSpout2, HIGH);
    digitalWrite(WaterSpout2_copy, HIGH);
    Serial.println(command == 'w' ? "ACK SMALL_REWARD" : "ACK LARGE_REWARD");
  } else if (command == 'p') {
    if (behaviorAir.active) { Serial.println("ERR AIR_BUSY"); return; }
    behaviorAir.active = true;
    behaviorAir.started = millis();
    behaviorAir.duration = SmallPunishSize;
    ++ManualNegativeNum;
    digitalWrite(Airpuff, HIGH);
    digitalWrite(Airpuff_copy, HIGH);
    Serial.println("ACK PUNISHMENT");
  } else if (command == 't') {
    behaviorBeginTone();
    Serial.println("ACK TONE");
  } else if (command == 'b' || command == 'r') {
    if (command == 'b') { behaviorBlue.active = false; digitalWrite(ShutterBlue, LOW); }
    else { behaviorRed.active = false; digitalWrite(ShutterRed, LOW); }
    Serial.println(command == 'b' ? "ACK BLUE_OPEN" : "ACK RED_OPEN");
  } else if (command == 'f') {
    if (behaviorStartPattern(behaviorBlue, ShutterBlue, BluePatternPulseMs, BluePatternPeriodMs, BluePatternPulseCount))
      Serial.println("ACK BLUE_PATTERN");
  } else if (command == 'F') {
    if (behaviorStartPattern(behaviorRed, ShutterRed, RedPatternPulseMs, RedPatternPeriodMs, RedPatternPulseCount))
      Serial.println("ACK RED_PATTERN");
  } else if (command == 'c') {
    if (behaviorWater.active || CalibrationRepeats == 0 || UnitRewardSize == 0) {
      Serial.println("ERR CALIBRATION_BUSY_OR_INVALID");
      return;
    }
    behaviorCalibrating = true;
    behaviorCalibrationDelivered = 0; // reset only when explicitly starting a new calibration
    behaviorCalibrationWaterOn = false;
    behaviorCalibrationChanged = millis();
    Serial.println("ACK CALIBRATION");
  }
}

void behaviorReadSerial() {
  // Bound work per loop so output timers still run under continuous input.
  for (byte n = 0; n < 16 && Serial.available() > 0; ++n) {
    behaviorHandleCommand((char)Serial.read());
  }
}
#endif
