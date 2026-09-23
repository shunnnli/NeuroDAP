// Bonsai to Arduino pulse stim
// Shun 03/27/2025

// ----- Parameters -----
double pulse_freq = 50.0;        // Hz
double pulse_duration = 5.0;     // ms (pulse width)
double pulse_interval;           // ms (delay between pulses)
unsigned long last_pulse_time = 0;

const byte triggerPin = 4;       // Input from Bonsai (digital HIGH to trigger)
const byte stimOutPin = 24;       // Output pin for stimulation pulses

bool isStimulating = false;

// ----- Sync Parameters -----
const byte Sync = 2;                  // Output pin for sync pulses
unsigned long TimerSync = 0;
unsigned long SyncPulseInterval = 100 + random(1, 400);
byte SyncNow = 0;

void setup() {
  pinMode(triggerPin, INPUT);    // Setup trigger pin as input
  pinMode(stimOutPin, OUTPUT);   // Setup stimulation pin as output
  digitalWrite(stimOutPin, HIGH); // Initialize output low
  pulse_interval = (1000.0 / pulse_freq) - pulse_duration;

  pinMode(Sync, OUTPUT);
  digitalWrite(Sync, LOW);
  TimerSync = millis();
}

void loop() {
  sync();

  // Check if the trigger line is HIGH
  if (digitalRead(triggerPin) == HIGH) {
    // If trigger is high, start patterned stimulation
    unsigned long currentTime = millis();

    if (currentTime - last_pulse_time >= pulse_interval + pulse_duration) {
      // Send pulse
      digitalWrite(stimOutPin, LOW);
      delay(pulse_duration);
      digitalWrite(stimOutPin, HIGH);

      last_pulse_time = millis(); // Update pulse time after turning off
    }
  } else {
    // If not stimulating, ensure output is low
    digitalWrite(stimOutPin, HIGH);
    last_pulse_time = millis(); // Reset timer to avoid burst after signal goes HIGH again
  }
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