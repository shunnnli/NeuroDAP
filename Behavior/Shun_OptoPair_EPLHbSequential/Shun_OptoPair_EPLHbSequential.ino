// Shun Li, 2022/11/10
// 2023/02/10
// 1. Add a 100ms delay between tone offset and airpuff onset
// 2. Reward trials no longer counter as trials
// 3. Tidied up the code more (remove block structure)

//20231016
// 1. Separate blue vs red opto delivery
// 2. Add code for paAIP2 stim
// 3. Change reward to pavlovian with 1000 delay

//20240502
// 1. Add option to choose between using blue or red to stimulate

//20240915
// Adapted to EP as cue and LHb stim as outcome
// Current try: replace airpuff with LHb stim

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
unsigned long SmallRewardSize = 5 * UnitRewardSize;
unsigned long BigRewardSize = 5 * UnitRewardSize;
unsigned long SmallPunishSize = 50;
unsigned long BigPunishSize = 200;

// Outcome probability params
// CHANGE BACK!!
int PairProbRange[2] = {1,40}; //40%
int StimOnlyProbRange[2] = {41,90}; //50%
int ToneOnlyProbRange[2] = {91,100}; //10%

int PunishProb = 90;
int OmissionProb = -1;
boolean OmitToneOnly = true;
boolean OmitStimOnly = false;

// Opto stim params
boolean RedStim = true; // 1 means using red to stim, 0 means using blue to stim
int TestAfter = 2000000; // Number of punish trials finished to enter test phase
int nTestPulse = 1; // Number of stim only trials in the test phase
int StimTotalPulseNum = 25; // number of pulses per pattern
unsigned long StimPulseDuration = 5; // Total duration of each pulse within a stimulation
unsigned long StimPulseFreq = 50; //For single 500ms pulse
unsigned long StimTotalDuration = 500; // Total duration of each stimulation
unsigned long StimPulseInterval = 0;

// Red opto stim params
// if RedStim if true, these will be rewritten to match opto stim params
int RedTotalPulseNum = 25; // number of pulses per pattern
unsigned long RedPulseDuration = 5; // Total duration of each pulse within a stimulation
unsigned long RedPulseFreq = 50; //For single 500ms pulse
unsigned long RedStimDuration = 500; // Total duration of each stimulation
unsigned long RedPulseInterval = 0;

// Blue opto stim params (basically a trigger pulse to matlab galvo.m)
// if RedStim if false, these will be rewritten to match opto stim params
int BlueTotalPulseNum = 15; // number of pulses per pattern
unsigned long BluePulseDuration = 10; // Total duration of each pulse within a stimulation
unsigned long BluePulseFreq = 30; //For single 500ms pulse
unsigned long BlueStimDuration = 500; // Total duration of each stimulation
unsigned long BluePulseInterval = 0;
// paAIP2 params
int paAIP2MaxTrial = 0;//10000;
unsigned long paAIP2StimInterval = 5000;
unsigned long paAIP2TimerPulse = 1000;

// Trial settings params
boolean pavlovian = false; // If false, omit punishment if animal blinks between CS-US
int minLicks = 2;

// Tone params
int LeftCueFreq = 3000;
int RightCueFreq = 12000;
unsigned long ToneDelayTime = 0; // start tone 250ms after

// Time params
unsigned long DelayTime = 0; // delay period between cue and outcome
unsigned long ReactionTime = 2000; // maximum reaction period (after cue) in ms
unsigned long TimeOutDuration = 10000; // time out duration in ms
unsigned long LongToneDuration = 1000; // Duration of the speaker tone in ms
unsigned long ShortToneDuration = 500;
unsigned long ToneDuration = 0; // Tone duration of each trial
unsigned long ITI1 = 2000;
unsigned long ITI2 = 4000;
unsigned long ITIMax = 3000000;
unsigned long ITIMin = 15000;
unsigned long ITIGracePeriod = 1000;
unsigned long ITI = 0; // ITI = random(ITI1,ITI2)

// Optotag pattern parameters (not used)
unsigned long ITIlaser = 25; // Time of stim in (ms)after cue
unsigned long ITIblue = 1000; // ITI of blue laser in ms
unsigned long PulseDurationBlue = 20; // duration of blue laser in ms
unsigned long PulseDurationRed = 2;
unsigned long ITIred = 10000; // duration of red laser in ms

// Block structure params (not used)
//int BlockLength = 10000; // Number of trials for each block
//boolean SwitchOnce = 0; // If 1, swtich after trials > SwitchAfter and only switch once; If 0, switch after BlockLength is reached
//int SwitchAfter = 35; // Number of trials before switching valence

//********** Params Initializtion ***********//
// Outcome related params
int trialRewardProb = 0;
int trialPunishProb = 0;
int trialOmissionProb = 0;
int trialFreeRewardProb = 0;
int trialRandomProb = 0;
boolean trialSecondCue = false;

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
const byte EyeClosure = 28; //Input from Bonsai (true if eye closes)
const byte ShutterBlue = 22; //1=blue shutter open, 0=closed
const byte ShutterRed = 24; //1=red shutter open, 0=closed


// Initialize real time variables //
char SerialInput = '0'; //for incoming serial data

// Constantly occuring stuffs
static int state = 0 ; // MAIN behavior state variable for running behavior task
int TrialType = 0; //0=LEFT, 1=RIGHT
unsigned long TimerSync = 0; //timer for non-periodic sync pulse
int SyncPulseInterval = 1000; //interval for non-periodic sync pulse
int SyncNow = 0; //current sync signal status

// Red Opto stim parameters
int RedPulseNum = 0;
unsigned long RedTimerPulse = 0;
int RedOptoNow = 0;
unsigned long RedOptoInterval = 0;

// Blue Opto stim parameters
int BluePulseNum = 0;
unsigned long BlueTimerPulse = 0;
int BlueOptoNow = 0;
unsigned long BlueOptoInterval = 0;

//Control shuttersound
unsigned long TimerShutterSound = 0;
unsigned long ShutterSoundInterval = 0;
int ShutterSoundNow = 0;

// Misc
boolean giveTest = false;
int remainingTest = 0;

// Trial structure related
int TrialNum = 0; //current trial number
// int BlockNum = 1; //current block number
// int TrialInBlock = 0; //current number of trials of this block
int PositiveNum = 0; //current positive outcome number
int NegativeNum = 0; //current negative outcome number
int ObtainedReward = 0; //reward follow by licking

// Detection related
int Lick = 0;
//int Blink = 0; //Whether there is eye closure between CS-US
int AnticipatoryLick = 0;
unsigned long LastLick = 0; //timestamp of last lick
//int LickLeft = 0; //left spout lick status
//int LickRight = 0; //right spout lick status
unsigned long Lick_Duration; //duration of current lick in ms
//int BlinkCount = 0;
int LickCount = 0;

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
int StimToneNum = 0;
int StimOnlyNum = 0;
int ToneOnlyNum = 0;
int OmissionNum = 0;

// boolean for printTrials()
int getReward = 0; //0: no reward; 1: small reward; 2: large reward
int getPunish = 0; //0: no punish; 1: small punish; 2: large punish
int getFreeReward = 0; //0: no free reward; 1: small free reward; 2: large free reward
int getLHbStim = 0;
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
  pinMode(EyeClosure,INPUT);
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
  // Initialize opto params
  RedTimerPulse = 0;
  RedOptoNow = 0;
  BlueTimerPulse = 0;
  BlueOptoNow = 0;

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
  Serial.println("Manual check: 1 -> reward; 2 -> punishment; 3 -> blue; 4 -> red");
  Serial.println("Laser shutter: 5 -> blue stim; 6 -> red stim");
  Serial.println("Water calibration: 7");
  Serial.println("Trial start/stop: 8 -> start; 9 -> end");
  Serial.println("-----------------------------------------------------------------");

  if (RedStim){
    // Red opto stim params
    RedTotalPulseNum = StimTotalPulseNum; // number of pulses per pattern
    RedPulseDuration = StimPulseDuration; // Total duration of each pulse within a stimulation
    RedPulseFreq = StimPulseFreq; //For single 500ms pulse
    RedStimDuration = StimTotalDuration; // Total duration of each stimulation
    RedPulseInterval = StimPulseInterval;

    // Print checks
    Serial.println("Stim color: red");
    Serial.print("Stim total pulse num: ");
    Serial.println(RedTotalPulseNum);
    Serial.print("Stim pulse duration: ");
    Serial.println(RedPulseDuration);
    Serial.print("Stim pulse freq: ");
    Serial.println(RedPulseFreq);
    Serial.print("Stim total duration: ");
    Serial.println(RedStimDuration);
    
  }else{
    // Blue opto stim params
    BlueTotalPulseNum = StimTotalPulseNum; // number of pulses per pattern
    BluePulseDuration = StimPulseDuration; // Total duration of each pulse within a stimulation
    BluePulseFreq = StimPulseFreq; //For single 500ms pulse
    BlueStimDuration = StimTotalDuration; // Total duration of each stimulation
    BluePulseInterval = StimPulseInterval;

    // Print checks
    Serial.println("Stim color: blue");
    Serial.print("Stim total pulse num: ");
    Serial.println(BlueTotalPulseNum);
    Serial.print("Stim pulse duration: ");
    Serial.println(BluePulseDuration);
    Serial.print("Stim pulse freq: ");
    Serial.println(BluePulseFreq);
    Serial.print("Stim total duration: ");
    Serial.println(BlueStimDuration);
  }
}


void loop() {
  sync(); //Non period sync pulse (1s width) generation
  lickDetection();
  //blinkDetection();
  opto();
  //randomShutterSound(ShutterBlue);

  // Check: if using blue to stim, no paAIP2 stim anymore
  if (!RedStim && paAIP2MaxTrial > 0){
    Serial.println("ERROR: should not use blue for stimulation if paAIP2 is on!");
    delay(1000000000000000);
  }

  // Give blue stim for paAIP2 before paAIP2MaxTrial
  if (millis() - paAIP2TimerPulse >= paAIP2StimInterval && TrialNum < paAIP2MaxTrial){
      paAIP2TimerPulse = millis();
      giveBlueOpto();
  }

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
      noTone(Speaker);
      ITI_start = millis();
      ITI_firstStart = millis();
      ITI = random(ITI1, ITI2);
      trialSecondCue = false;
      //updateTrials(); //If update trial here, ITI following tone belongs to next trial
      printTrials();
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
        updateTrials();
        LickCount = 0;
        if (giveTest) {
          trialPunishProb = 1;
          if (OmitStimOnly){
            trialOmissionProb = -100; // no outcome  
          }else{
            trialOmissionProb = 100;
          }
          trialFreeRewardProb = 0;

          if (RedStim){giveRedOpto();}
          else{giveBlueOpto();}
          ToneDuration = ShortToneDuration;
          StimOnlyNum += 1;

        } else {
          trialPunishProb = random(1,101);
          trialOmissionProb = random(1,101);
          trialRandomProb = random(1,101); // trial type (pair, stimOnly, toneOnly)

          if (trialPunishProb > PunishProb) {
            ToneDuration = 0;
          } else {
            ToneDuration = ShortToneDuration;
            if (trialRandomProb >= PairProbRange[0] && trialRandomProb <= PairProbRange[1]){
              if (RedStim){giveRedOpto();}
              else{giveBlueOpto();}
              trialSecondCue = true;
              StimToneNum += 1;
            }else if (trialRandomProb >= StimOnlyProbRange[0] && trialRandomProb <= StimOnlyProbRange[1]){
              if (RedStim){giveRedOpto();}
              else{giveBlueOpto();}
              ToneDuration = ShortToneDuration;
              StimOnlyNum += 1;
              if (OmitStimOnly){
                trialOmissionProb = -100; // no outcome for tone only trials  
              }
            }else if (trialRandomProb >= ToneOnlyProbRange[0] && trialRandomProb <= ToneOnlyProbRange[1]){
              if (OmitToneOnly){
                trialOmissionProb = -100; // no outcome for tone only trials  
              }
              trialSecondCue = true;
              ToneOnlyNum += 1;
            }
          }
        }

        Cue_start = millis();
        TrialType = 0;
        getFreeReward = 0;
        printTrials();
        state = 3;
      }

      if (Lick == 1 && ENL) {
        ITI_start = millis();
      }
      break;

    //state 3; second cue on
    case SecondCueOn:
      if (trialSecondCue && millis() - Cue_start >= ToneDelayTime) {
        Cue_start = millis();
        tone(Speaker, LeftCueFreq);
        digitalWrite(SpeakerLeft_copy, HIGH);
        state = 4;
      } else if (!trialSecondCue && millis() - Cue_start >= ToneDelayTime) {
        Cue_start = millis();
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
        if ((Cue_off - Cue_start) > ToneDuration){
          state = 5; 
        }
      }
      break;

    //state 5: Turn on solenoid
    case SolenoidOn:
      if (trialPunishProb > PunishProb) {
        getReward = 2;
        getPunish = 0;
        getLHbStim = 0;
        giveReward();
        Punish_start = 0;
      } else {
        if (!pavlovian && LickCount >= minLicks){
          Miss += 1;
          if (trialOmissionProb > OmissionProb){
            getReward = 0;
            getPunish = 0;
            getLHbStim = 1;
            Reward_start = 0;
            Punish_start = 0;
            state = 6;
          } else {
            getReward = 0;
            getPunish = 0;
            getLHbStim = 0;
            state = 6;
          }
        } else if (millis()-Cue_start > ReactionTime){
          if (pavlovian){
            if (trialOmissionProb > OmissionProb){
              getReward = 0;
              getPunish = 0;
              getLHbStim = 1;
              Reward_start = 0;
              Punish_start = 0;
              state = 6;
            } else {
              getReward = 0;
              getPunish = 0;
              getLHbStim = 0;
              state = 6;
            }
          } else {
            Hit += 1;
            if (trialOmissionProb > OmissionProb){
              getReward = 0;
              getPunish = 0;
              getLHbStim = 1;
              Reward_start = 0;
              Punish_start = 0;
              state = 6;
            } else {
              getReward = 0;
              getPunish = 0;
              getLHbStim = 0;
              state = 6;
            }
          }
        }
      }

      // Record reward/punishment delivery for current trial
      trialReward = getReward;
      trialPunish = getPunish;
      trialFreeReward = getFreeReward;
      break;

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

      if (getLHbStim == 1){
        giveBlueOpto();
        getLHbStim = 0;
      }

      if (getReward == 0 && getPunish == 0 && getLHbStim == 0) {
        printTrials();
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
        printTrials(); // print after time out ends
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
      Serial.println(millis());
      Serial.print("Total reward: ");
      Serial.println(PositiveNum);
      Serial.print("Total punishment: ");
      Serial.println(NegativeNum);
      Serial.print("Total stim only trials: ");
      Serial.println(StimOnlyNum);
      Serial.print("Total tone only trials: ");
      Serial.println(ToneOnlyNum);
      Serial.print("Total stim&tone trials: ");
      Serial.println(StimToneNum);
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

    //Serial.print("I received: ");
    //Serial.println(SerialInput);

    if (SerialInput == '1' && LeftOutcomeButton == 0) { // dispense left reward
      Serial.println("Entered 1: reward");
      //tone(Speaker, LeftCueFreq);
      //digitalWrite(SpeakerLeft_copy, HIGH);
      //digitalWrite(ShutterRed, HIGH);
      //delay(ToneDuration);
      //noTone(Speaker);
      //digitalWrite(SpeakerLeft_copy, LOW);
      //digitalWrite(ShutterRed, LOW);
      LeftOutcomeButton = 1;
      digitalWrite(WaterSpout2, HIGH);
      digitalWrite(WaterSpout2_copy, HIGH);
      OutcomeSize = BigRewardSize;
      LeftOutcomeTimer = millis();
    }

    if (SerialInput == '2' && RightOutcomeButton == 0) {
      Serial.println("Entered 2: punishment");
      tone(Speaker, LeftCueFreq);
      digitalWrite(SpeakerLeft_copy, HIGH);
      delay(500);
      noTone(Speaker);
      digitalWrite(SpeakerLeft_copy, LOW);
      RightOutcomeButton = 1;
      digitalWrite(Airpuff, HIGH);
      digitalWrite(Airpuff_copy, HIGH);
      OutcomeSize = SmallPunishSize;
      RightOutcomeTimer = millis();
    }

    if (SerialInput == '3') {
      Serial.println("Entered 3: Deliver blue stim");
      // giveBlueOpto();
      digitalWrite(ShutterBlue,LOW);
    }

    if (SerialInput == '4') {
      Serial.println("Entered 4: Deliver red stim");
//      giveRedOpto();
      digitalWrite(ShutterRed,LOW);
    }

    if (SerialInput == '5') { //blue pulsing pattern
      Serial.println("Entered 5: blue stim");
//      for (int i = 0; i < 20; i++) {
//        Serial.print("Red stim #");
//        Serial.println(i + 1);
//        for (int j = 0; j < 5; j++) {
//          digitalWrite(ShutterRed, HIGH);
//          delay(PulseDurationRed);
//          //delay(4);
//          digitalWrite(ShutterRed, LOW);
//          delay(50 - PulseDurationRed);
//          //delay(46);
//        }
//        delay(9500);
//        //delay(19000); // for j == 20, 20Hz
//        //delay(18000); // for j == 40, 20Hz
//        //delay(18000); // for j == 40, 20Hz
//      }
      giveBlueOpto();
      Serial.println("5: blue pulsing finished");
    }

    if (SerialInput == '6') { //red pulsing pattern
      Serial.println("Entered 6: red pulsing");
//      for (int i = 0; i < 50; i++) {
//        digitalWrite(ShutterBlue, HIGH);
//        delay(PulseDurationBlue);
//        digitalWrite(ShutterBlue, LOW);
//        delay(ITIblue - PulseDurationBlue);
//      }
      giveRedOpto();
      Serial.println("6: red pulsing finished");
    }

    if (SerialInput == '7'){
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
      SyncPulseInterval = 100 + random(400); // random sync pulse interval between 1~2s
    } else {
      TimerSync = millis();
      digitalWrite(Sync, HIGH);
      SyncNow = 1;
      SyncPulseInterval = 50;
    }
  }
}

//********************************************************************************************//
// Assigned upcoming red opto delivery
void giveRedOpto(){
  if (RedTotalPulseNum != 1){
    RedPulseNum = (RedStimDuration/1000.0) * RedPulseFreq;
    RedPulseInterval = (1000.0/RedPulseFreq) - RedPulseDuration;
  } else{
    RedPulseNum = RedTotalPulseNum;
    RedPulseInterval = 5;
  }
  
  if (RedPulseInterval <= 0 && RedPulseNum > 1){
    RedPulseInterval = 5;
    Serial.println("Negative RedPulseInterval: reset to 5ms");
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
  if (RedPulseNum > 0 && millis() - RedTimerPulse >= RedOptoInterval){
    if (RedOptoNow == 1){
      RedTimerPulse = millis();
      digitalWrite(ShutterRed, HIGH);
      RedOptoNow = 0;
      RedPulseNum -= 1;
      RedOptoInterval = RedPulseInterval;
    } else {
      RedTimerPulse = millis();
      digitalWrite(ShutterRed, LOW);
      RedOptoNow = 1;
      RedOptoInterval = RedPulseDuration;
    }
  }

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
// Deliver random shutter sound
void randomShutterSound(int ShutterColor) {
  if (millis() - TimerShutterSound >= ShutterSoundInterval) {
    if (ShutterSoundNow == 1) {
      TimerShutterSound = millis();
      digitalWrite(ShutterColor, LOW);
      ShutterSoundNow = 0;
      ShutterSoundInterval = 100 + random(1, 900); // random opto interval between 10~20s
      //Serial.print("Next ShutterSound after ");
      //Serial.print(ShutterSoundInterval / 1000.0);
      //Serial.println("s");
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
//void blinkDetection(){
//  if (digitalRead(EyeClosure) == 1) {
//    if (Blink == 0) {
//      Blink = 1;
//      BlinkCount += 1;
//      Serial.print("Trial: ");
//      Serial.print(TrialNum);
//      Serial.print("\t");
//
//      Serial.print("Blink Detected");
//      Serial.print("\t");
//
//      Serial.print("Time: ");
//      Serial.println(millis() / 1000.0);
//      
//    }
//  }
//  if (digitalRead(EyeClosure) == 0) {
//    if (Blink == 1) {
//      Blink = 0;
//    }
//  }
//}

//********************************************************************************************//
void updateTrials() {
  // Update trial and block
  if (!giveTest) {
    TrialNum = TrialNum + 1;

    if (remainingTest == 0 && (TrialNum - 1) % TestAfter == 0) {
      giveTest = true;
      remainingTest = nTestPulse;
      Serial.println("New test started");
    }else if (trialReward > 0){
      // Only punishment trials are counted as trials (reward does not count)
      TrialNum -= 1;
    }

  } else {
    //Serial.println(remainingTest);
    // Give test stim if needed
    remainingTest -= 1;
    if (remainingTest == 0) {
      giveTest = false;
      Serial.println("Stim test finished");
    }
  }
}

//********************************************************************************************//
void printTrials() {
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

    if (trialRandomProb >= PairProbRange[0] && trialRandomProb <= PairProbRange[1]) {
      Serial.print("Cue start (Pair #");
      Serial.print(StimToneNum);
      Serial.print(")");
      Serial.print("\t");
    }else if (trialRandomProb >= StimOnlyProbRange[0] && trialRandomProb <= StimOnlyProbRange[1]){
      Serial.print("Cue start (Stim only #");
      Serial.print(StimOnlyNum);
      Serial.print(")");
      Serial.print("\t");
    }else if (trialRandomProb >= ToneOnlyProbRange[0] && trialRandomProb <= ToneOnlyProbRange[1]){
      Serial.print("Cue start (Tone only #");
      Serial.print(ToneOnlyNum);
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
      Serial.print("Small reward");
      Serial.print("\t");
    } else if (trialReward == 2) {
      Serial.print("Big reward");
      Serial.print("\t");
    } else if (trialPunish == 1) {
      if (giveTest) {
        Serial.print("Small punish (stim test #");
        Serial.print(nTestPulse - remainingTest + 1);
        Serial.print(")");
        Serial.print("\t");
      } else {
        Serial.print("Small punish");
        Serial.print("\t");
      }
    } else if (trialPunish == 2) {
      if (giveTest) {
        Serial.print("Big punish (stim test #");
        Serial.print(nTestPulse - remainingTest + 1);
        Serial.print(")");
        Serial.print("\t");
      } else {
        Serial.print("Big punish");
        Serial.print("\t");
      }
    } else if (trialReward == 0 && trialPunish == 0) {
      if (giveTest) {
        Serial.print("Omission (stim test #");
        Serial.print(nTestPulse - remainingTest + 1);
        Serial.print(")");
        Serial.print("\t");
      } else {
        Serial.print("Omission");
        Serial.print("\t");
      }
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
