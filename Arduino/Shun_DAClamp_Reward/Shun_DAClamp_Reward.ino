//Shun Li, 10/28/2022
// 20230210
// 1. Separate from FlexibleLearning_OneTone
// 2. Remove block structure
// 3. Tidied up the code

//20230722
// 1. Add anticipatory licks
// 2. If anticipatory licks + choice lick >= 2, get reward

//20231016
// 1. Separate blue vs red opto delivery
// 2. Add code for paAIP2 stim
// 3. Change reward to pavlovian with 1000 delay]

//20240502
// 1. Add option to choose between using blue or red to stimulate

//20260429
// 1. Remove all laser/opto/paAIP2 code
// 2. Randomly choose LeftCueFreq (50%) vs RightCueFreq (50%) per trial
// 3. Separate OmissionProb for Left and Right cue

#define Idle 0
#define ITI_State 1
#define FirstCueOn 2
#define SecondCueOn 3
#define FirstCueOff 4
#define SolenoidOn 5
#define SolenoidOff 6
#define TurnOffLickLeft 7
#define RestartClock 8
#define TimeOut 9

#include <math.h>

//********** User settings ***********//
// Set up parameters for the behavior
boolean ENL = true; // whether ITI is ENL
unsigned long UnitRewardSize = 20; // reward size of 1ul
unsigned long SmallRewardSize = 3 * UnitRewardSize;
unsigned long BigRewardSize = 8 * UnitRewardSize;
unsigned long SmallPunishSize = 50;
unsigned long BigPunishSize = 200;

// Cue selection (50/50 left vs right)
int LeftCueProb = 100; // probability (1-100) of using LeftCueFreq

// Outcome probability params (per cue)
int LeftOmissionProb = 10;  // omission prob for Left cue trials
int RightOmissionProb = 10; // omission prob for Right cue trials
int FreeRewardProb = 100;

// Trial settings params
boolean pavlovian = true; // If true, make reward pavlovian; if false, reward is operant
int minLicks = 2; // min amount of licks within response window to get an reward
int minLicks_pav = 4; // min amount of licks to get big reward for pavlovian task

// Tone params
int LeftCueFreq = 11000;
int RightCueFreq = 6000;
unsigned long ShortToneDuration = 500;
unsigned long ToneDelayTime = 400; // start tone 250ms after for clamping

// Time params
unsigned long DelayTime = 1500; // delay period between cue and outcome
unsigned long ReactionTime = 1500; // maximum reaction period (after cue) in ms
unsigned long TimeOutDuration = 10000; // time out duration in ms
unsigned long ITI1 = 2000;
unsigned long ITI2 = 4000;
unsigned long ITIMax = 3000000;
unsigned long ITIMin = 15000;
unsigned long ITIGracePeriod = 1000;
unsigned long ITI = 0; // ITI = random(ITI1,ITI2)

//********** Params Initializtion ***********//
// Outcome related params
int trialOmissionProb = 0;
int trialFreeRewardProb = 0;
int trialCueProb = 0;
boolean trialUseLeftCue = true; // which cue this trial uses

// Trial counters
int TrialNum = 0;        // current trial number
int LeftCueNum = 0;      // total left-cue trials
int RightCueNum = 0;     // total right-cue trials
int PositiveNum = 0;     // current positive outcome number
int NegativeNum = 0;     // current negative outcome number

// Misc
static int state = 0 ; // MAIN behavior state variable for running behavior task
unsigned long ToneDuration = 0; // Tone duration of each trial

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


// Initialize real time variables //
char SerialInput = '0'; //for incoming serial data

// Constantly occuring stuffs
unsigned long TimerSync = 0; //timer for non-periodic sync pulse
int SyncPulseInterval = 1000; //interval for non-periodic sync pulse
int SyncNow = 0; //current sync signal status

// Detection related
int Lick = 0;
unsigned long LastLick = 0; //timestamp of last lick
unsigned long Lick_Duration; //duration of current lick in ms
int LickCount = 0;

// Outcome related
unsigned long OutcomeSize = 20; // current trial's port reward size
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

//Summary response variables
int Hit = 0;
int Miss = 0;
int FalseAlarm = 0;
int CorrectReject = 0;
int GoNum = 0;
int NoGoNum = 0;

// boolean for printTrials()
int getReward = 0; //0: no reward; 1: small reward; 2: large reward
int getPunish = 0; //0: no punish; 1: small punish; 2: large punish
int getFreeReward = 0; //0: no free reward; 1: small free reward; 2: large free reward
int trialReward = getReward;
int trialPunish = getPunish;
int trialFreeReward = getFreeReward;

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

  Serial.println("-----------------------------------------------------------------");
  Serial.println("Manual check: 1 -> reward; 2 -> tone; 3 -> blue shutter; 4 -> red shutter");
  Serial.println("Water calibration: 7");
  Serial.println("Trial start/stop: 8 -> start; 9 -> end");
  Serial.print("LeftCueFreq: ");
  Serial.println(LeftCueFreq);
  Serial.print("RightCueFreq: ");
  Serial.println(RightCueFreq);
  Serial.print("LeftCueProb: ");
  Serial.println(LeftCueProb);
  Serial.print("LeftOmissionProb: ");
  Serial.println(LeftOmissionProb);
  Serial.print("RightOmissionProb: ");
  Serial.println(RightOmissionProb);
  Serial.println("-----------------------------------------------------------------");
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
      digitalWrite(ShutterBlue, HIGH);
      digitalWrite(ShutterRed, HIGH);
      noTone(Speaker);
      ITI_start = millis();
      ITI_firstStart = millis();
      ITI = random(ITI1, ITI2);
      printTrials(state, trialReward, trialPunish);
      state = 2;
      break;

    //state 2: Turn on Cue
    case FirstCueOn:
      Cue_start = millis();
      Actual_ITI = Cue_start - ITI_firstStart;
      Current_ITI = Cue_start - ITI_start;
      trialITIMin = random(ITIMin, ITIMin + ITIGracePeriod);
      trialITIMax = random(ITIMax - ITIGracePeriod, ITIMax);

      if (Current_ITI > ITI && Actual_ITI > trialITIMin) {
        TrialNum += 1;
        LickCount = 0;

        // Choose left vs right cue (50/50 by LeftCueProb)
        trialCueProb = random(1, 101);
        if (trialCueProb <= LeftCueProb) {
          trialUseLeftCue = true;
          LeftCueNum += 1;
        } else {
          trialUseLeftCue = false;
          RightCueNum += 1;
        }

        trialOmissionProb = random(1, 101);
        trialFreeRewardProb = random(1, 101);

        ToneDuration = ShortToneDuration;
        Cue_start = millis();
        getFreeReward = 0;
        printTrials(state, trialReward, trialPunish);
        state = 3;
      }

      if (Lick == 1 && ENL) {
        ITI_start = millis();
      }
      break;

    // state 3: Second cue (tone) on
    case SecondCueOn:
      if (millis() - Cue_start >= ToneDelayTime) {
        Cue_start = millis();
        if (trialUseLeftCue) {
          tone(Speaker, LeftCueFreq);
          digitalWrite(SpeakerLeft_copy, HIGH);
        } else {
          tone(Speaker, RightCueFreq);
          digitalWrite(SpeakerRight_copy, HIGH);
        }
        state = 4;
      }
      break;

    //state 4: Turn off cue
    case FirstCueOff:
      Cue_off = millis();
      if ((Cue_off - Cue_start) > ToneDuration) {
        noTone(Speaker);
        digitalWrite(SpeakerLeft_copy, LOW);
        digitalWrite(SpeakerRight_copy, LOW);
        if (millis() - Cue_start > ToneDuration) {
          state = 5;
        }
      }
      break;

    //state 5: Turn on solenoid after delay period (1 sec)
    case SolenoidOn: {
      // Pick the cue-specific omission threshold for this trial
      int currentOmissionProb = trialUseLeftCue ? LeftOmissionProb : RightOmissionProb;

      if (!pavlovian && LickCount >= minLicks) {
        Hit += 1;
        if (trialOmissionProb > currentOmissionProb) {
          getReward = 2;
          getPunish = 0;
          giveReward();
          Punish_start = 0; // so it automatically shuts the solenoid down
        } else {
          getReward = 0;
          getPunish = 0;
          state = 6;
        }
      } else if (millis() - Cue_start > ReactionTime) {
        if (pavlovian) {
          if (LickCount >= minLicks_pav) {
            if (trialOmissionProb > currentOmissionProb) {
              getReward = 2;
              getPunish = 0;
              giveReward();
              Punish_start = 0;
            } else {
              getReward = 0;
              getPunish = 0;
              state = 6;
            }
          } else {
            if (trialOmissionProb > currentOmissionProb) {
              getReward = 1;
              getPunish = 0;
              giveReward();
              Punish_start = 0;
            } else {
              getReward = 0;
              getPunish = 0;
              state = 6;
            }
          }
        } else {
          Miss += 1;
          if (trialFreeRewardProb <= FreeRewardProb) {
            if (trialOmissionProb > currentOmissionProb) {
              getFreeReward = 1;
              getReward = 1;
              getPunish = 0;
              giveReward();
              Punish_start = 0;
            } else {
              getReward = 0;
              getPunish = 0;
              state = 6;
            }
          } else {
            state = 9; //Time out
            printTrials(state, trialReward, trialPunish);
          }
        }
      }

      // Record reward/punishment delivery for current trial
      trialReward = getReward;
      trialPunish = getPunish;
      trialFreeReward = getFreeReward;
      break;
    }

    //state 6: Turn off solenoid
    case SolenoidOff:
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
        printTrials(state, trialReward, trialPunish);
        if (trialFreeReward > 0) {
          state = 9; //go to time out trial after free reward
        } else {
          state = 1; //go to state 1 if all reward and punish are given
        }
      }
      break;

    //state 9: time out -> return to ITI
    case TimeOut:
      Timeout_start = millis();
      if ((Timeout_start - Cue_off - DelayTime) > TimeOutDuration) {
        printTrials(state, trialReward, trialPunish); // print after time out ends
        state = 1;
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
      Serial.println(millis() / 1000.0);
      Serial.print("Total reward: ");
      Serial.println(PositiveNum);
      Serial.print("Total punishment: ");
      Serial.println(NegativeNum);
      Serial.print("Total left cue trials: ");
      Serial.println(LeftCueNum);
      Serial.print("Total right cue trials: ");
      Serial.println(RightCueNum);
      Serial.print("Go: ");
      Serial.println(GoNum);
      Serial.print("NoGoNum: ");
      Serial.println(NoGoNum);
      Serial.print("Hit: ");
      Serial.println(Hit);
      Serial.print("Miss: ");
      Serial.println(Miss);
      Serial.print("FalseAlarm: ");
      Serial.println(FalseAlarm);
      Serial.print("CorrectReject: ");
      Serial.println(CorrectReject);
      state = 0;
      noTone(Speaker);
      End = millis();
    }
  }

  if (Serial.available() > 0) {
    // read the incoming byte:
    SerialInput = Serial.read();

    if (SerialInput == '1' && LeftOutcomeButton == 0) { // dispense left reward
      Serial.println("Entered 1: reward");
      LeftOutcomeButton = 1;
      digitalWrite(WaterSpout2, HIGH);
      digitalWrite(WaterSpout2_copy, HIGH);
      OutcomeSize = SmallRewardSize;
      LeftOutcomeTimer = millis();
    }

    if (SerialInput == '2') {
      Serial.println("Entered 2: Tone");
      tone(Speaker, LeftCueFreq);
      digitalWrite(SpeakerLeft_copy, HIGH);
      delay(ShortToneDuration);
      noTone(Speaker);
      digitalWrite(SpeakerLeft_copy, LOW);
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
void giveReward() {
  Reward_start = millis();
  digitalWrite(WaterSpout2, HIGH);
  digitalWrite(WaterSpout2_copy, HIGH);
  On = millis();
  PositiveNum += 1;
  state = 6;
}

//********************************************************************************************//
void givePunishment() {
  Punish_start = millis();
  digitalWrite(Airpuff, HIGH);
  digitalWrite(Airpuff_copy, HIGH);
  On = millis();
  NegativeNum += 1;
  state = 6;
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
      Lick = 0;
      Serial.print("Trial: ");
      Serial.print(TrialNum);
      Serial.print("\t");

      Serial.print("Lick Detected");
      Serial.print("\t");

      Lick_Duration = millis() - Lick_Start;
      Serial.print("Time: ");
      Serial.println(millis() / 1000.0);
      LastLick = millis();
      LickCount += 1;
    }
  }
  // Lick Detection END //
}

//********************************************************************************************//
void printTrials(int state, int trialReward, int trialPunish) {
  // ITI state
  if (state == 1) {
    Serial.print("Trial: ");
    Serial.print(TrialNum);
    Serial.print("\t");

    Serial.print("ITI = ");
    Serial.print(ITI / 1000.0);
    Serial.print("\t");

    Serial.print("Time: ");
    Serial.println(ITI_start / 1000.0);
  }
  else if (state == 2) {
    Serial.print("Trial: ");
    Serial.print(TrialNum);
    Serial.print("\t");

    if (trialUseLeftCue) {
      Serial.print("Cue start (Left #");
      Serial.print(LeftCueNum);
      Serial.print(")");
      Serial.print("\t");
    } else {
      Serial.print("Cue start (Right #");
      Serial.print(RightCueNum);
      Serial.print(")");
      Serial.print("\t");
    }

    Serial.print("Time: ");
    Serial.print(Cue_start / 1000.0);
    Serial.println("\t");
  }
  else if (state == 5) {
    Serial.print("Trial: ");
    Serial.print(TrialNum);
    Serial.print("\t");
    Serial.println("I don't think state 5 have printTrials()");

  }
  else if (state == 6) {
    Serial.print("Trial: ");
    Serial.print(TrialNum);
    Serial.print("\t");

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
    } else if (trialReward == 0 && trialPunish == 0) {
      Serial.print("Reward block: omission");
      Serial.print("\t");
    }

    Serial.print("Time: ");
    Serial.println(millis() / 1000.0);
  }
  else if (state == 9) {
    Serial.print("Trial: ");
    Serial.print(TrialNum);
    Serial.print("\t");

    Serial.print("Time out");
    Serial.print("\t");

    Serial.print("Time: ");
    Serial.println(millis() / 1000.0);
  }
}
