// Important Parameters for the task that can be changed ////

unsigned long RewardSizeLeft = 10;
unsigned long RewardSizeRight = 25;

int num_repeat = 200; // number of free rewards per command (push button or serial command)


const byte Speaker = 47; //speaker output pin
const byte LickDetect1 = 18; //left lick detection
const byte LickDetect2 = 19; //right lick detection
const byte WaterSpout = 8; //left spout solenoid
const byte WaterSpout2 = 9; //right spout solenoid
const byte WaterSpout_copy = 6; //copy left spout solenoid for data recording device
const byte WaterSpout2_copy = 7; //copy right spout solenoid for data recording device
const byte Button = 53; //start-stop button
const byte Button1 = 11; //left spout manual activation
const byte Button2 = 12; //right spout manual activation
const byte CameraTrigger = 21; //camera trigger
const byte StartStop = 30; //start-stop session output pin
const byte Laser = 51; //start-stop session output pin
const byte ShutterBlue = 22; //1=blue shutter open, 0=closed
const byte ShutterRed = 24; //1=red shutter open, 0=closed

int ButtonStatus1 = 0;
int ButtonStatus2 = 0;
unsigned long DummyButton1 = 0;
unsigned long DummyButton2 = 0;

int GoCueFreq = 12000;
int NoGoCueFreq = 3000;

char SerialInput = '0'; //for incoming serial data

// GUI manual-control settings. This hardware-check sketch has no airpuff or task.
#define Idle 0
static int state = Idle;
unsigned long UnitRewardSize = 25;
unsigned long SmallRewardSize = 10;
unsigned long BigRewardSize = 25;
unsigned long ShortToneDuration = 500;
#define LeftCueFreq GoCueFreq
#define CalibrationRepeats num_repeat
unsigned long CalibrationIntervalMs = 100;
unsigned long BluePatternPulseMs = 5;
unsigned long BluePatternPeriodMs = 50;
unsigned int BluePatternPulseCount = 25;
unsigned long RedPatternPulseMs = 5;
unsigned long RedPatternPeriodMs = 50;
unsigned int RedPatternPulseCount = 25;
int LeftOutcomeButton = 0;
int RightOutcomeButton = 0;
int ManualPositiveNum = 0;
int ManualNegativeNum = 0;
unsigned int checkToneSteps = 0;
unsigned long checkToneChanged = 0;
#define BEHAVIOR_HAS_AIR 0
#define BEHAVIOR_HAS_TASK 0
#define BEHAVIOR_HAS_TONE_COPIES 0
#define BEHAVIOR_SHUTTER_OPEN HIGH
#define BEHAVIOR_SHUTTER_CLOSED LOW
#define BEHAVIOR_LEGACY_REWARD_BIG 1
#define BEHAVIOR_LEGACY_SECOND_TONE 0
#define BEHAVIOR_CUSTOM_LEGACY
bool behaviorHandleLegacy(char command);
#include "BehaviorManual.h"

void setup()
{
  Serial.begin(115200);

  pinMode(Button, INPUT);

  pinMode(WaterSpout, OUTPUT);
  pinMode(WaterSpout2, OUTPUT);
  pinMode(WaterSpout_copy, OUTPUT);
  pinMode(WaterSpout2_copy, OUTPUT);
  pinMode(Button1, INPUT);
  pinMode(Button2, INPUT);
  pinMode(Speaker, OUTPUT);
  noTone(Speaker);
  pinMode(LickDetect1, INPUT);
  pinMode(LickDetect2, INPUT);
  pinMode(ShutterBlue, OUTPUT);
  pinMode(ShutterRed, OUTPUT);

  digitalWrite(WaterSpout, LOW);
  digitalWrite(WaterSpout2, LOW);
  digitalWrite(WaterSpout_copy, LOW);
  digitalWrite(WaterSpout2_copy, LOW);
  digitalWrite(ShutterBlue, LOW);
  digitalWrite(ShutterRed, LOW);
  noTone(Speaker);

  behaviorCapabilities();
  behaviorReportState();

  Serial.println("type 1 for left reward / 2 for right reward");
  Serial.println("for laser shutter: type 3 for blue / 4 for red");
}

// Legacy 1/2 and physical buttons still calibrate the left/right spout.
void checkCalibration(bool left, unsigned long interval) {
  if (behaviorCalibrating || behaviorWater.active || checkToneSteps) {
    Serial.println("ERR BUSY stop the calibration/manual stimulus first");
    return;
  }
  behaviorHandleCommand('c');
  if (behaviorCalibrating) {
    behaviorCalibrationValve = left ? WaterSpout : WaterSpout2;
    behaviorCalibrationCopy = left ? WaterSpout_copy : WaterSpout2_copy;
    behaviorCalibrationSize = left ? RewardSizeLeft : RewardSizeRight;
    behaviorCalibrationInterval = interval;
  }
}

bool behaviorHandleLegacy(char command) {
  if (command == '1' || command == '2') {
    checkCalibration(command == '1', 100);
    return true;
  }
  if (command == '5') {
    if (behaviorCalibrating || checkToneSteps) {
      Serial.println("ERR BUSY stop calibration/tone first");
      return true;
    }
    checkToneSteps = 100; // Original 50 alternating pairs at 100 ms per tone.
    checkToneChanged = millis();
    tone(Speaker, GoCueFreq);
    Serial.println("ACK TONE");
    return true;
  }
  // Keep GUI stimuli from competing with the legacy tone sequence.
  if (checkToneSteps && command != 'x' && command != '9') {
    Serial.println("ERR BUSY tone sequence in progress");
    return true;
  }
  return false;
}

void loop() {
  behaviorReadSerial();
  behaviorUpdateManual();
  if (checkToneSteps && millis() - checkToneChanged >= 100) {
    checkToneChanged = millis();
    if (--checkToneSteps) tone(Speaker, checkToneSteps % 2 ? NoGoCueFreq : GoCueFreq);
    else noTone(Speaker);
  }
  int left = digitalRead(Button1);
  int right = digitalRead(Button2);
  if (left && !ButtonStatus1) checkCalibration(true, 50);
  if (right && !ButtonStatus2) checkCalibration(false, 100);
  ButtonStatus1 = left;
  ButtonStatus2 = right;
  if (digitalRead(LickDetect1) == LOW) Serial.println("Left lick detected");
  if (digitalRead(LickDetect2) == LOW) Serial.println("Right lick detected");
}

void behaviorStartTask() {} // No automated task in this hardware-check sketch.
void behaviorStopTask() {
  state = Idle;
  checkToneSteps = 0;
  behaviorStopOutputs();
}
