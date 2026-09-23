// Shun Li, 2022/11/10
// Random outcome task: give water reward, punishment, or tone randomly

//20230210: tidied up code, renamed to RandomOutcome

//20260625
// 1. Remove extra outcome code
// 2. Match sync pulse timing to Shun_DAClamp_Reward

#define Idle 0
#define ITI_State 1
#define OutcomeDelivery 2
#define EndOutcomeDelivery 3
//#define SecondCueOn 4
#define SolenoidOn 5
#define SolenoidOff 6
#define TurnOffLickLeft 7
#define RestartClock 8
#define TimeOut 9
#include <math.h>

//********** User settings ***********//
// Set up for the behavior
boolean ENL = true; // whether ITI is ENL
unsigned long UnitRewardSize = 15; // reward size of 1ul
unsigned long SmallRewardSize = 3 * UnitRewardSize;
unsigned long BigRewardSize = 8 * UnitRewardSize;
unsigned long SmallPunishSize = 100;
unsigned long BigPunishSize = 200;

// Outcome probability params
int RewardProbRange[2] = {0, 50}; //{0, 40};
int PunishProbRange[2] = {51, 90}; //{41, 50};
int ToneProbRange[2] = {91, 100}; //{51, 60};

// Time dependent params
unsigned long ShortToneDuration = 500;
unsigned long LongToneDuration = 1000;
unsigned long eventDelayTime = 500; // delay between ITI satisfaction and event delivery
unsigned long ITI1 = 2000;
unsigned long ITI2 = 4000;
unsigned long ITIMax = 200000; //200000;
unsigned long ITIMin = 15000; //15000;
unsigned long ITIGracePeriod = 1000;
unsigned long ITI = 0;

// Tone params
int LeftCueFreq = 3000;
int RightCueFreq = 12000;

//********** Params Initializtion ***********//
// Input output pin description //
const byte Sync = 2; // non-periodic sync pulse
const byte Speaker = 47; //speaker output pin
const byte SpeakerLeft_copy = 48; //speaker output step copy - Left tone
const byte SpeakerRight_copy = 49; //speaker output step copy - Right tone
const byte LickDetectLeft = 18; //left lick detection
const byte LickDetectRight = 19; //right lick detection
const byte WaterSpout = 8; //left spout solenoid
const byte WaterSpout2 = 9; //right spout solenoid
const byte WaterSpout_copy = 6; //copy left spout solenoid for data recording device
const byte WaterSpout2_copy = 7; //copy right spout solenoid for data recording device
const byte Airpuff = 32; //airpuff valve
const byte Airpuff_copy = 34; //airpuff valve copy for data receiving device
const byte ShutterBlue = 22; //1=blue shutter closed, 0=open
const byte ShutterRed = 24; //1=red shutter closed, 0=open

// Misc
char SerialInput = '0'; //for incoming serial data
unsigned long ToneDuration = 0;
int trialRandomProb = 0;

// Constantly occuring stuffs
static int state = 0 ; // MAIN behavior state variable for running behavior task
unsigned long TimerSync = 0; //timer for non-periodic sync pulse
int SyncPulseInterval = 1000; //interval for non-periodic sync pulse
int SyncNow = 0; //current sync signal status
//randomReward()
unsigned long TimerReward = 0;
unsigned long RewardInterval = 0;
int RewardNow = 0;
//randomPunish()
unsigned long TimerPunish = 0;
unsigned long PunishInterval = 0;
int PunishNow = 0;

// Trial structure related
int PositiveNum = 0; //current positive outcome number
int NegativeNum = 0;
int ManualPositiveNum = 0;
int ManualNegativeNum = 0;
int toneNum = 0;
int TrialNum = 0;

// boolean for printTrials()
int getReward = 0; //0: no reward; 1: small reward; 2: large reward
int getPunish = 0; //0: no punish; 1: small punish; 2: large punish
int getFreeReward = 0; //0: no free reward; 1: small free reward; 2: large free reward
int trialReward = getReward;
int trialPunish = getPunish;
int trialFreeReward = getFreeReward;

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

// Timestamp related
unsigned long Lick_Start = 0; //timestamp for current lick
unsigned long ITI_start = 0; //timestamp for beginning of ITI
unsigned long Current_ITI = 0; // current ITI (reset by licks)
unsigned long Actual_ITI = 0; // actual elapsed time from last trial (not reset by licks)
unsigned long ITI_firstStart = 0;
unsigned long EventDelay_start = 0; //timestamp for event delay after ITI is satisfied
unsigned long trialITIMin = 0;
unsigned long trialITIMax = 0;
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

void setup()
{
  Serial.begin(115200);

  pinMode(Sync, OUTPUT);
  pinMode(WaterSpout, OUTPUT);
  pinMode(WaterSpout2, OUTPUT);
  pinMode(Airpuff, OUTPUT);
  pinMode(WaterSpout_copy, OUTPUT);
  pinMode(WaterSpout2_copy, OUTPUT);
  pinMode(Airpuff_copy, OUTPUT);
  pinMode(LickDetectLeft, INPUT);
  pinMode(LickDetectRight, INPUT);
  pinMode(Speaker, OUTPUT);
  pinMode(SpeakerLeft_copy, OUTPUT);
  pinMode(SpeakerRight_copy, OUTPUT);
  pinMode(ShutterBlue, OUTPUT);
  pinMode(ShutterRed, OUTPUT);

  // Initialize sync params
  Start = millis();
  TimerSync = millis();
  state = 0;
  SyncNow = 0;
  // Initialize random reward params
  TimerReward = 0;
  RewardNow = 1;
  // Initialize random punish params
  TimerPunish = 0;
  PunishNow = 1;

  digitalWrite(Sync, LOW);
  digitalWrite(WaterSpout, LOW);
  digitalWrite(WaterSpout2, LOW);
  digitalWrite(Airpuff, LOW);
  digitalWrite(WaterSpout_copy, LOW);
  digitalWrite(WaterSpout2_copy, LOW);
  digitalWrite(Airpuff_copy, LOW);
  noTone(Speaker);
  digitalWrite(ShutterBlue, HIGH);
  digitalWrite(ShutterRed, HIGH);
  randomSeed(analogRead(3));

  Serial.println("---------------------------------RandomOutcome--------------------------------");
  Serial.println("Manual check: 1 -> reward; 2 -> punishment; 3 -> blue shutter; 4 -> red shutter");
  Serial.println("Water calibration: 7");
  Serial.println("Trial start/stop: 8 -> start; 9 -> end");
  Serial.println("---------------------------------RandomOutcome--------------------------------");
}


void loop() {
  sync(); //Non period sync pulse (1s width) generation
  lickDetection();

  switch (state) {
    //state 0: Idle state until Start button pushed
    case Idle:
      if (SerialInput == '8') {
        Start = millis();
        End = 0;
        Serial.print("TASK STARTED AT ");
        Serial.print("\t");
        Serial.println(millis());
        state = 1;
        digitalWrite(WaterSpout, LOW);
        digitalWrite(WaterSpout2, LOW);
        digitalWrite(Airpuff, LOW);
        digitalWrite(WaterSpout_copy, LOW);
        digitalWrite(WaterSpout2_copy, LOW);
        digitalWrite(Airpuff_copy, LOW);
        digitalWrite(ShutterBlue, HIGH);
        digitalWrite(ShutterRed, HIGH);
        noTone(Speaker);
      }
      break;

    //state 1: Determine the intertrial interval
    case ITI_State:
      ITI_start = millis();
      ITI_firstStart = millis();
      EventDelay_start = 0;
      ITI = random(ITI1, ITI2);
      state = 2;
      Serial.print("ITI: ");
      Serial.println(ITI / 1000.0);
      break;


    //state 2: select what to deliver
    case OutcomeDelivery:
      Now = millis();
      Actual_ITI = Now - ITI_firstStart;
      Current_ITI = Now - ITI_start;
      trialITIMin = random(ITIMin - ITIGracePeriod, ITIMin + ITIGracePeriod);
      trialITIMax = random(ITIMax - ITIGracePeriod, ITIMax + ITIGracePeriod);

      if (Lick == 1 && ENL) {
        ITI_start = millis();
        EventDelay_start = 0;
        break;
      }

      if ((Current_ITI > ITI && Actual_ITI > trialITIMin) || EventDelay_start > 0) {
        if (EventDelay_start == 0) {
          EventDelay_start = millis();
          TrialNum += 1;

          Serial.print("Trial: ");
          Serial.print(TrialNum);
          Serial.print("\t");
          Serial.print("ITI finished");
          Serial.print("\t");
          Serial.print("Time: ");
          Serial.println(EventDelay_start / 1000.0);
        }

        if (millis() - EventDelay_start >= eventDelayTime) {
          Cue_start = millis();
          trialRandomProb = random(101);
          //Serial.println(RewardProbRange[1]);
          if (trialRandomProb >= RewardProbRange[0] && trialRandomProb <= RewardProbRange[1]) {
            Reward_start = millis();
            digitalWrite(WaterSpout2, HIGH);
            digitalWrite(WaterSpout2_copy, HIGH);
            PositiveNum += 1;
            Punish_start = 0;
            getReward = 2;

            Serial.print("Reward: ");
            Serial.print(PositiveNum);
            Serial.print("\t");
            Serial.print("Time: ");
            Serial.println(millis() / 1000.0);

          } else if (trialRandomProb >= PunishProbRange[0] && trialRandomProb <= PunishProbRange[1]) {
            Punish_start = millis();
            digitalWrite(Airpuff, HIGH);
            digitalWrite(Airpuff_copy, HIGH);
            NegativeNum += 1;
            Reward_start = 0;
            getPunish = 2;

            Serial.print("Punish: ");
            Serial.print(NegativeNum);
            Serial.print("\t");
            Serial.print("Time: ");
            Serial.println(millis() / 1000.0);

          } else if (trialRandomProb >= ToneProbRange[0] && trialRandomProb <= ToneProbRange[1]) {
            tone(Speaker, LeftCueFreq);
            digitalWrite(SpeakerLeft_copy, HIGH);
            delay(ShortToneDuration);
            noTone(Speaker);
            digitalWrite(SpeakerLeft_copy, LOW);
            toneNum += 1;

            Serial.print("Tone: ");
            Serial.print(toneNum);
            Serial.print("\t");
            Serial.print("Time: ");
            Serial.println(millis() / 1000.0);
          }
          EventDelay_start = 0;
          state = 3;
        }
      } else {
        EventDelay_start = 0;
      }

      // Record reward/punishment delivery for current trial
      trialReward = getReward;
      trialPunish = getPunish;
      trialFreeReward = getFreeReward;
      break;

    //state 3: end reward delivery
    case EndOutcomeDelivery:
      Outcome_off = millis();
      if (getReward == 1) {
        if ((Outcome_off - Reward_start) >= SmallRewardSize) {
          digitalWrite(WaterSpout, LOW);
          digitalWrite(WaterSpout_copy, LOW);
          digitalWrite(WaterSpout2, LOW);
          digitalWrite(WaterSpout2_copy, LOW);
          getReward = 0;
        }
      } else if (getReward == 2) {
        if ((Outcome_off - Reward_start) >= BigRewardSize) {
          digitalWrite(WaterSpout, LOW);
          digitalWrite(WaterSpout_copy, LOW);
          digitalWrite(WaterSpout2, LOW);
          digitalWrite(WaterSpout2_copy, LOW);
          getReward = 0;
        }
      }

      if (getPunish == 1) {
        if ((Outcome_off - Punish_start) >= SmallPunishSize) {
          digitalWrite(Airpuff, LOW);
          digitalWrite(Airpuff_copy, LOW);
          getPunish = 0;
        }
      } else if (getPunish == 2) {
        if ((Outcome_off - Punish_start) >= BigPunishSize) {
          digitalWrite(Airpuff, LOW);
          digitalWrite(Airpuff_copy, LOW);
          getPunish = 0;
        }
      }

      if (getReward == 0 && getPunish == 0) {
        printTrials();
        state = 1; //go to state 1 if all reward and punish are given
      }
      break;

    default:
      state = 1;
      break;
  }
  // END OF SWITCH STRCUTURE //


  // Ending Task //
  if (SerialInput == '9') {
    if (End == 0) {
      Serial.print("TASK ENDED AT ");
      Serial.print("\t");
      Serial.println(millis());
      Serial.print("Total reward: ");
      Serial.println(PositiveNum);
      Serial.print("Total punish: ");
      Serial.println(NegativeNum);
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
      ManualPositiveNum += 1;
      Serial.print("Entered 1: manual reward #");
      Serial.println(ManualPositiveNum);
      LeftOutcomeButton = 1;
      digitalWrite(WaterSpout2, HIGH);
      digitalWrite(WaterSpout2_copy, HIGH);
      OutcomeSize = BigRewardSize;
      LeftOutcomeTimer = millis();
    }

    if (SerialInput == '2' && RightOutcomeButton == 0) {
      Serial.print("Entered 2: manual punishment #");
      Serial.println(ManualNegativeNum);
      RightOutcomeButton = 1;
      digitalWrite(Airpuff, HIGH);
      digitalWrite(Airpuff_copy, HIGH);
      OutcomeSize = SmallPunishSize;
      RightOutcomeTimer = millis();
    }
    if (SerialInput == '3') {
      Serial.println("Entered 3: Open blue shutter");
      digitalWrite(ShutterBlue, LOW);
    }
    if (SerialInput == '4') {
      Serial.println("Entered 4: Open red shutter");
      digitalWrite(ShutterRed, LOW);
    }
    if (SerialInput == '7') {
      Serial.println("Deliver 200 reward for calibration");
      int num_repeat = 200;
      for (int i = 0; i < num_repeat; i++) {
        delay(100);
        digitalWrite(WaterSpout2, HIGH);
        digitalWrite(WaterSpout2_copy, HIGH);
        delay(UnitRewardSize);
        digitalWrite(WaterSpout2, LOW);
        digitalWrite(WaterSpout2_copy, LOW);
      }
      Serial.print("Finished delivery x");
      Serial.println(num_repeat);
    }

  }

  //Left Button outcome
  if ((millis() - LeftOutcomeTimer > OutcomeSize) && LeftOutcomeButton == 1) {
    digitalWrite(WaterSpout2, LOW);
    digitalWrite(WaterSpout2_copy, LOW);
    Serial.print("Manual reward");
    Serial.print("\t");
    Serial.print(LeftOutcomeTimer / 1000.0);
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
    Serial.print(RightOutcomeTimer / 1000.0);
    Serial.print("\t");
    Serial.println(OutcomeSize);
    RightOutcomeButton = 0;
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

//********************************************************************************************//
void randomReward() {
  if (millis() - TimerReward >= RewardInterval) {
    if (RewardNow == 1) {
      TimerReward = millis();
      digitalWrite(WaterSpout2, LOW);
      digitalWrite(WaterSpout2_copy, LOW);
      RewardNow = 0;
      RewardInterval = 5000 + random(10000); // random reward interval between 5 - 15s

      Serial.print("Next reward after ");
      Serial.print(RewardInterval / 1000.0);
      Serial.println("s");
    } else {
      TimerReward = millis();
      digitalWrite(WaterSpout2, HIGH);
      digitalWrite(WaterSpout2_copy, HIGH);
      RewardNow = 1;
      PositiveNum += 1;
      RewardInterval = BigRewardSize;
      printTrials();
    }
  }
}

//********************************************************************************************//
void randomPunish() {
  if (millis() - TimerPunish >= PunishInterval) {
    if (PunishNow == 1) {
      TimerPunish = millis();
      digitalWrite(Airpuff, LOW);
      digitalWrite(Airpuff_copy, LOW);
      PunishNow = 0;
      PunishInterval = 20000 + random(1, 30000); // random reward interval between 20~50s
      Serial.print("Next punishment after ");
      Serial.print(PunishInterval / 1000.0);
      Serial.println("s");
    } else {
      TimerPunish = millis();
      digitalWrite(Airpuff, HIGH);
      digitalWrite(Airpuff_copy, HIGH);
      PunishNow = 1;
      NegativeNum += 1;
      PunishInterval = SmallPunishSize;
      printTrials();
    }
  }
}

//********************************************************************************************//
void lickDetection() {
  // Lick Detection //
  if (digitalRead(LickDetectRight) == 0) {
    if (Lick == 0) {
      Lick = 1;
      Lick_Start = millis();
    }
  }
  if (digitalRead(LickDetectRight) == 1) {
    if (Lick == 1) {
      delay(1);
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
  // ITI state
  if (state == 3) {
    if (trialReward == 1) {
      Serial.print("Reward block: small reward");
      Serial.print("\t");
    } else if (trialReward == 2) {
      Serial.print("Reward block: big reward");
      Serial.print("\t");
    } else if (trialPunish == 1) {
      Serial.print("Reward block: small punish");
      Serial.print("\t");
    } else if (trialPunish == 2) {
      Serial.print("Reward block: big punish");
      Serial.print("\t");
    }
    Serial.print("Time: ");
    Serial.println(millis() / 1000.0);
  }
}
