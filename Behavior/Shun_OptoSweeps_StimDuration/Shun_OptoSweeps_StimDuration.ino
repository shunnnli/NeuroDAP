// Shun_OptoSweeps_Frequency
// Shun Li, 2023/07/18
// Modified from code by Ricardo López, 2022/12/15

// Give pulses of opto stimulation at fixed define pulse duration and frequency
// but different frequency

#define Idle 0 //Defining constants and assigning values. These are variables with fixed values defined before the program runs.
#define ITI_State 1
#define OptoSweeps 2
//#define SweepInterval 3
//#define SecondCueOn 4
//#define SolenoidOn 5
//#define SolenoidOff 6
//#define TurnOffLickLeft 7
//#define RestartClock 8
//#define TimeOut 9

#include <math.h>
// Define opto stim parameters
int PatternNum = 3; //Total number of different frequencies
unsigned long StimDurationRange[3] = {5,15,25}; // in seconds
unsigned long BluePulseFreq = 50; // in Hz
unsigned long BluePulseDuration = 5; //Total duration of each pulse
// Initialize opto stim params
unsigned long BlueStimDuration = 0; //Current pulse frequency in the current sweep
unsigned long BluePulseInterval = 0;
unsigned long BlueTotalPulseNum = 0;
// Blue Opto stim parameters
int BluePulseNum = 0;
unsigned long BlueTimerPulse = 0;
int BlueOptoNow = 0;
unsigned long BlueOptoInterval = 0;
// Initialize pattern counting
unsigned long PatternCount[3] = {0,0,0};


// Set up parameters for the behavior
unsigned long SmallRewardSize = 2 * 8; //Unsigned long varaibles are extended size variables for number storage and store 32 bits
unsigned long BigRewardSize = 5 * 8; //Here we assign reward and punishment sizes
unsigned long SmallPunishSize = 50;
unsigned long BigPunishSize = 100;

// Time params
unsigned long ITI1 = 40000;
unsigned long ITI2 = 40000;
unsigned long ITI = 0; // ITI = random(ITI1,ITI2)

// Input output pin description //Constants prevent specific object/method()/variable to modify data
const byte Sync = 2; // non-periodic sync pulse
const byte Speaker = 47; //speaker output pin
const byte SpeakerLeft_copy = 48; //speaker output step copy - Left tone
const byte SpeakerRight_copy = 49; //speaker output step copy - Right tone
const byte LickDetect1 = 18; //left lick detection
const byte LickDetect2 = 19; //right lick detection
const byte WaterSpout = 9; //left spout solenoid
const byte WaterSpout2 = 8; //right spout solenoid
const byte WaterSpout_copy = 6; //copy left spout solenoid for data recording device
const byte WaterSpout2_copy = 7; //copy right spout solenoid for data recording device
const byte Airpuff = 32; //airpuff valve
const byte Airpuff_copy = 34; //airpuff valve copy for data receiving device
const byte ShutterBlue = 22; //1=blue shutter open, 0=closed
const byte ShutterRed = 24; //1=red shutter open, 0=closed [This may be important for laser stimulation]

// Initialize real time variables //
char SerialInput = '0'; //for incoming serial data

// Constantly occuring stuffs
static int state = 0 ; // MAIN behavior state variable for running behavior task
unsigned long TimerSync = 0; //timer for non-periodic sync pulse
int SyncPulseInterval = 1000; //interval for non-periodic sync pulse
int SyncNow = 0; //current sync signal status
// Opto stim parameters
unsigned long TimerPulse = 0;
int OptoNow = 0;
unsigned long OptoInterval = 0;
boolean OptoDelivering = false;
// Random shutter sound params
unsigned long TimerShutterSound = 0;
unsigned long ShutterSoundInterval = 0;
int ShutterSoundNow = 0;

// Detection related
int Lick = 0;
unsigned long LastLick = 0; //timestamp of last lick
unsigned long Lick_Duration; //duration of current lick in ms

// Outcome related
unsigned long OutcomeSize = 20; // current trial's port reward size (updated by OutcomeSizeLeft or OutcomeSizeRight)
int LeftOutcomeButton = 0; //left reward button status
int RightOutcomeButton = 0; //right reward button status
unsigned long LeftOutcomeTimer = 0; //timer for left reward button
unsigned long RightOutcomeTimer = 0; //timer for right reward button
int trialRewardProb = 0;

// Timestamp related
unsigned long Lick_Start = 0; //timestamp for current lick
unsigned long ITI_start = 0; //timestamp for beginning of ITI
unsigned long Cue_start = 0; //timestamp for cue onset
unsigned long Cue_off = 0; //timestamp for cue off
unsigned long Reward_start = 0; //timestamp for start of reward
unsigned long Punish_start = 0; //timestamp for start of punishment
unsigned long Outcome_off = 0; //timestamp for solenoid off
unsigned long Timeout_start = 0; //timestamp for timeout
unsigned long On; //timestamp for reward delivery (solenoid on)
unsigned long Start = 0; //timestamp for starting the session (used for triggering camera)
unsigned long Now = 0; //current timesatmp (used for triggering camera)
unsigned long End = 0; //timestamp for ending the session

void setup() { //Setup function called when sketch starts. The setup function is used to initialize variables, pin modes, etc.

  Serial.begin(115200); //Set the data rate in bits per second. 115200 is default

  pinMode(Sync, OUTPUT); //Configure specified pins as inputs or outputs
  pinMode(WaterSpout, OUTPUT);
  pinMode(WaterSpout2, OUTPUT);
  pinMode(Airpuff, OUTPUT);
  pinMode(WaterSpout_copy, OUTPUT);
  pinMode(WaterSpout2_copy, OUTPUT);
  pinMode(Airpuff_copy, OUTPUT);
  pinMode(LickDetect1, INPUT);
  pinMode(LickDetect2, INPUT);
  pinMode(Speaker, OUTPUT);
  pinMode(SpeakerLeft_copy, OUTPUT);
  pinMode(SpeakerRight_copy, OUTPUT);
  pinMode(ShutterBlue, OUTPUT);
  pinMode(ShutterRed, OUTPUT);

  // Initialize sync params
  Start = millis(); //Number of milliseconds passed since program starts
  TimerSync = millis();
  state = 0; //Sync parameters are defined above beginning at L62
  SyncNow = 0;
  // Initialize opto params
  BlueTimerPulse = 0;
  BlueOptoNow = 0;

  digitalWrite(Sync, LOW); //Set pins off
  digitalWrite(WaterSpout, LOW);
  digitalWrite(WaterSpout2, LOW);
  digitalWrite(Airpuff, LOW);
  digitalWrite(WaterSpout_copy, LOW);
  digitalWrite(WaterSpout2_copy, LOW);
  digitalWrite(Airpuff_copy, LOW);
  noTone(Speaker);
  digitalWrite(ShutterBlue, HIGH);
  digitalWrite(ShutterRed, HIGH);

  Serial.println("-----------------------------------------------------------------");
  Serial.println("Manual check: 1 -> reward; 2 -> punishment; 3 -> blue; 4 -> red");
  Serial.println("Laser shutter: 5 -> red stim; 6 -> blue stim");
  Serial.println("Trial start/stop: 8 -> start; 9 -> end");
  Serial.println("-----------------------------------------------------------------");
}


void loop() {
  sync(); //Non period sync pulse (1s width) generation
  lickDetection();
  opto(); 

  switch (state) {
    //state 0: Idle state until Start button pushed
    case Idle: //Recall Idle = 0 so we begin in our first state here
      if (SerialInput == '8') {
        Start = millis();
        End = 0;
        Serial.print("TASK STARTED AT ");
        Serial.print("\t");
        Serial.println(millis());
        digitalWrite(WaterSpout, LOW);
        digitalWrite(WaterSpout2, LOW);
        digitalWrite(Airpuff, LOW);
        digitalWrite(WaterSpout_copy, LOW);
        digitalWrite(WaterSpout2_copy, LOW);
        digitalWrite(Airpuff_copy, LOW);
        noTone(Speaker);
        digitalWrite(ShutterBlue, HIGH);
        state = 2;
        ITI_start = -3000;
      }
      break;

    //state 1: Determine the intertrial interval
    case ITI_State:
      digitalWrite(WaterSpout, LOW);
      digitalWrite(WaterSpout2, LOW);
      digitalWrite(Airpuff, LOW);
      digitalWrite(WaterSpout_copy, LOW);
      digitalWrite(WaterSpout2_copy, LOW);
      digitalWrite(Airpuff_copy, LOW);
      noTone(Speaker);
      ITI_start = millis();
      ITI = random(ITI1,ITI2);
      //printTrials(state, trialReward, trialPunish);
      state = 2;
      break;

    //state 2: Determine the intertrial interval
    case OptoSweeps:
      if (millis()-ITI_start > ITI){
        // Randomly choose the stim duration
        int randChoice = random(0,3);
        BlueStimDuration = StimDurationRange[randChoice] * 1000;
        PatternCount[randChoice] += 1; 

        // Deliver stim
        giveBlueOpto();
        Serial.print("Finished: #");
        Serial.print(PatternCount[randChoice]);
        Serial.print(" opto sweep at ");
        Serial.print(BlueStimDuration/1000.0);
        Serial.println(" s");
        state = 1;
      }
      
      break;

    // state 3: interval between each sweep pattern
//    case SweepInterval:
//      if (PatternNum == CurrentPattern && CurrentRepeat == RepeatPerPattern){
//        Serial.println("Finished: all opto sweeps");
//        state = 1;
//      }
//      state = 1;
//      break;
  }
  // END OF SWITCH STRCUTURE //


  // Ending Task //
  if (SerialInput == '9') {
    if (End == 0) {
      Serial.print("TASK ENDED AT ");
      Serial.print("\t");
      Serial.println(millis());
      state = 0;
      noTone(Speaker);
      End = millis();
    }
  }

  if (Serial.available() > 0) {
    // read the incoming byte:
    SerialInput = Serial.read();

    //Serial.print("I received: ");
    //Serial.println(SerialInput);

    if (SerialInput == '1' && LeftOutcomeButton == 0) { // dispense left reward
      Serial.println("Entered 1: reward");
      LeftOutcomeButton = 1;
      digitalWrite(WaterSpout2, HIGH);
      digitalWrite(WaterSpout2_copy, HIGH);
      OutcomeSize = BigRewardSize;
      LeftOutcomeTimer = millis();
    }

    if (SerialInput == '2' && RightOutcomeButton == 0) {
      Serial.println("Entered 2: punishment");
      RightOutcomeButton = 1;
      digitalWrite(Airpuff, HIGH);
      digitalWrite(Airpuff_copy, HIGH);
      OutcomeSize = SmallPunishSize;
      RightOutcomeTimer = millis();
    }

    if (SerialInput == '3') {
      digitalWrite(ShutterBlue, HIGH);
      //delay(2000);
      //digitalWrite(ShutterBlue, LOW);
    }

    if (SerialInput == '4') {
      digitalWrite(ShutterRed, HIGH);
      //delay(2000);
      //digitalWrite(ShutterRed, LOW);
    }

    // if (SerialInput == '5') { //red pulsing pattern
    //   Serial.println("Entered 5: red pulsing");
    //   for (int i = 0; i < 20; i++) { //For loop runs 20 times
    //     Serial.print("Red stim #");
    //     Serial.println(i + 1); //Print red laser stimulus number
    //     for (int j = 0; j < 5; j++) { //For loop runs 5 times
    //       digitalWrite(ShutterRed, HIGH); //Red shutter opens
    //       delay(PulseDurationRed); //Shutter remains open for a delay of 2ms
    //       //delay(4);
    //       digitalWrite(ShutterRed, LOW); //Red shutter closes
    //       delay(50 - PulseDurationRed); //48ms delay
    //       //delay(46);
    //     }
    //     delay(9500);
    //     //delay(19000); // for j == 20, 20Hz
    //     //delay(18000); // for j == 40, 20Hz
    //     //delay(18000); // for j == 40, 20Hz
    //   }
    //   Serial.println("5: red pulsing finished");
    // }

    // if (SerialInput == '6') { //blue pulsing pattern
    //   Serial.println("Entered 6: blue pulsing");
    //   for (int i = 0; i < 50; i++) {
    //     digitalWrite(ShutterBlue, HIGH);
    //     delay(PulseDurationBlue);
    //     digitalWrite(ShutterBlue, LOW);
    //     delay(ITIblue - PulseDurationBlue);
    //   }
    //   Serial.println("6: blue pulsing finished");
    // }

  }

  //Left Button outcome
  if ((millis() - LeftOutcomeTimer > OutcomeSize) && LeftOutcomeButton == 1) {
    digitalWrite(WaterSpout2, LOW);
    digitalWrite(WaterSpout2_copy, LOW);
    Serial.print("Manual reward");
    Serial.print("\t");
    Serial.print(LeftOutcomeTimer);
    Serial.print("\t");
    Serial.println(OutcomeSize);
    LeftOutcomeButton = 0;
  }

  //Right button outcome
  if ((millis() - RightOutcomeTimer) > OutcomeSize && RightOutcomeButton == 1) {
    digitalWrite(Airpuff, LOW);
    digitalWrite(Airpuff_copy, LOW);
    Serial.print("Manual punishment");
    Serial.print("\t");
    Serial.print(RightOutcomeTimer);
    Serial.print("\t");
    Serial.println(OutcomeSize);
    RightOutcomeButton = 0;
  }
}

//********************************************************************************************//
void sync() { // Void function governing sync pulses?
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

//********************************************************************************************//
// Assigned upcoming blue opto delivery
void giveBlueOpto(){
  if (BlueTotalPulseNum != 1){
    BluePulseNum = (BlueStimDuration/1000.0) * BluePulseFreq;
    BluePulseInterval = (1000.0/BluePulseFreq) - BluePulseDuration;
  } else {
    BluePulseNum = BlueTotalPulseNum;
    BluePulseInterval = 5;
  }
  
  if (BluePulseInterval <= 0 && BluePulseNum > 1){
    BluePulseInterval = 5;
    Serial.println("Negative BluePulseInterval: reset to 5ms");
  }
}

//********************************************************************************************//
// Execute opto delivery 
void opto(){
//  if (RedPulseNum > 0 && millis() - RedTimerPulse >= RedOptoInterval){
//    if (RedOptoNow == 1){
//      RedTimerPulse = millis();
//      digitalWrite(ShutterRed, HIGH);
//      RedOptoNow = 0;
//      RedPulseNum -= 1;
//      RedOptoInterval = RedPulseInterval;
//    } else {
//      RedTimerPulse = millis();
//      digitalWrite(ShutterRed, LOW);
//      RedOptoNow = 1;
//      RedOptoInterval = RedPulseDuration;
//    }
//  }

  if (BluePulseNum > 0 && millis() - BlueTimerPulse >= BlueOptoInterval){
    if (BlueOptoNow == 1){
      BlueTimerPulse = millis();
      digitalWrite(ShutterBlue, HIGH);
      BlueOptoNow = 0;
      BluePulseNum -= 1;
      BlueOptoInterval = BluePulseInterval;
    } else {
      BlueTimerPulse = millis();
      digitalWrite(ShutterBlue, LOW);
      BlueOptoNow = 1;
      BlueOptoInterval = BluePulseDuration;
    }
  }
}

//********************************************************************************************//
void randomShutterSound(const byte ShutterColor) {
  if (millis() - TimerShutterSound >= ShutterSoundInterval) {
    if (ShutterSoundNow == 1) {
      TimerShutterSound = millis();
      digitalWrite(ShutterColor, LOW);
      ShutterSoundNow = 0;
      ShutterSoundInterval = 100 + random(1, 900); // random opto interval between 10~20s
    } else {
      TimerShutterSound = millis();
      digitalWrite(ShutterColor, HIGH);
      ShutterSoundNow = 1;
      ShutterSoundInterval = 500;
    }
  }
}

//********************************************************************************************//
void lickDetection() {
  // Lick Detection //
  if (digitalRead(LickDetect1) == 0 || digitalRead(LickDetect2) == 0) { //Or conditional
    if (Lick == 0) {
      Lick = 1;
      Lick_Start = millis(); //Assign starting time to lick
    }
  }
  if (digitalRead(LickDetect1) == 1 && digitalRead(LickDetect2) == 1) { //And conditional
    if (Lick == 1) {
      Lick = 0;

      Serial.print("Lick Detected");
      Serial.print("\t");

      Lick_Duration = millis() - Lick_Start;
      Serial.print("Time: ");
      Serial.println(millis() / 1000.0);
      LastLick = millis();
    }
  }
  // Lick Detection END //
}

//********************************************************************************************//
void printTrials() {
  
}
