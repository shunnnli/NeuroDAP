// Bonsai to Arduino pulse stim
// Shun 03/27/2025

// ----- Parameters -----
double pulse_freq = 50.0;        // Hz
double pulse_duration = 5.0;     // ms (pulse width)
unsigned long last_pulse_time = 0;

const byte triggerPin = 4;       // Input from Bonsai (digital HIGH to trigger)
const byte stimOutPin = 24;       // Output pin for stimulation pulses

bool isStimulating = false;

// ----- Sync Parameters -----
const byte Sync = 2;                  // Output pin for sync pulses
unsigned long TimerSync = 0;
unsigned long SyncPulseInterval = 100 + random(1, 400);
byte SyncNow = 0;

// GUI protocol: this sketch has a red stimulation output and an external-trigger task.
#define Idle 0
static int state = Idle; // Start task arms the Bonsai trigger after boot/reset.
#define ShutterRed stimOutPin
#define BEHAVIOR_HAS_WATER 0
#define BEHAVIOR_HAS_AIR 0
#define BEHAVIOR_HAS_TONE 0
#define BEHAVIOR_HAS_BLUE 0
#define BEHAVIOR_LEGACY_REWARD_BIG 1
#define BEHAVIOR_LEGACY_SECOND_TONE 0
#define RedPatternPulseMs ((unsigned long)pulse_duration)
#define RedPatternPeriodMs (pulse_freq > 0 ? (unsigned long)(1000.0 / pulse_freq) : 0UL)
unsigned int RedPatternPulseCount = 25;
#include "BehaviorManual.h"

void setup() {
  Serial.begin(115200);
  pinMode(triggerPin, INPUT);    // Setup trigger pin as input
  pinMode(stimOutPin, OUTPUT);   // Setup stimulation pin as output
  digitalWrite(stimOutPin, HIGH); // Active-low output starts closed

  pinMode(Sync, OUTPUT);
  digitalWrite(Sync, LOW);
  TimerSync = millis();
  behaviorCapabilities();
  behaviorReportState();
}

void loop() {
  behaviorReadSerial();
  behaviorUpdateManual();
  sync();
  if (state == Idle) return;

  // Nonblocking external-trigger pulses keep End task responsive during a pulse.
  unsigned long now = millis();
  if (digitalRead(triggerPin) != HIGH) {
    digitalWrite(stimOutPin, HIGH);
    isStimulating = false;
    last_pulse_time = now;
  } else if (isStimulating) {
    if (now - last_pulse_time >= pulse_duration) {
      digitalWrite(stimOutPin, HIGH);
      isStimulating = false;
      last_pulse_time = now;
    }
  } else if (pulse_freq > 0 && pulse_duration > 0 &&
             pulse_duration <= 1000.0 / pulse_freq &&
             now - last_pulse_time >= 1000.0 / pulse_freq) {
    digitalWrite(stimOutPin, LOW);
    isStimulating = true;
    last_pulse_time = now;
  }
}

void behaviorStartTask() {
  isStimulating = false;
  last_pulse_time = millis();
  state = 1;
}

void behaviorStopTask() {
  state = Idle;
  isStimulating = false;
  behaviorStopOutputs();
}

//********************************************************************************************//
void sync() {
  if (millis() - TimerSync >= SyncPulseInterval) {
    if (SyncNow == 1) {
      TimerSync = millis();
      digitalWrite(Sync, LOW);
      SyncNow = 0;
      SyncPulseInterval = 100 + random(1, 400); // random sync pulse interval between 1~2s
    } else {
      TimerSync = millis();
      digitalWrite(Sync, HIGH);
      SyncNow = 1;
      SyncPulseInterval = 50;
    }
  }
}
