// Included after the complete .ino; tests its real loop and manual-control code.
void sendByte(char command) { Serial.input.push_back(command); loop(); }
void advance(unsigned long milliseconds) { fakeTime += milliseconds; loop(); }
void assertClosed() {
  assert(pins[WaterSpout]==LOW && pins[WaterSpout2]==LOW && pins[Airpuff]==LOW);
  assert(pins[WaterSpout_copy]==LOW && pins[WaterSpout2_copy]==LOW && pins[Airpuff_copy]==LOW);
  assert(pins[ShutterBlue]==HIGH && pins[ShutterRed]==HIGH);
  assert(pins[SpeakerLeft_copy]==LOW && pins[SpeakerRight_copy]==LOW && !toneOn);
}
int main() {
  setup();
  assert(Serial.output.find("BEHAVIOR_CONTROLS 1")!=std::string::npos);
  // Standardized rewards use distinct durations in both sketches.
  sendByte('w'); assert(behaviorWater.duration==SmallRewardSize && pins[WaterSpout2]==HIGH);
  advance(SmallRewardSize); assert(pins[WaterSpout2]==LOW);
  sendByte('W'); assert(behaviorWater.duration==BigRewardSize);
  advance(BigRewardSize);
  sendByte('p'); assert(pins[Airpuff]==HIGH && !toneOn);
  advance(SmallPunishSize); assert(pins[Airpuff]==LOW);
  sendByte('t'); assert(toneOn);
  advance(ShortToneDuration); assert(!toneOn);
  // Legacy mappings remain sketch-specific.
  sendByte('1'); assert(behaviorWater.duration==(BEHAVIOR_LEGACY_REWARD_BIG ? BigRewardSize : SmallRewardSize));
  sendByte('x');
  sendByte('2');
  assert(BEHAVIOR_LEGACY_REWARD_BIG ? behaviorAir.active : toneOn);
  sendByte('x');
  sendByte('b'); assert(pins[ShutterBlue]==LOW);
  Serial.output.clear(); sendByte('?');
  assert(Serial.output.find("BEHAVIOR_STATE task=0 blue=1 red=0")!=std::string::npos);
  sendByte('B'); assert(pins[ShutterBlue]==HIGH);
  sendByte('r'); assert(pins[ShutterRed]==LOW);
  sendByte('R'); assert(pins[ShutterRed]==HIGH);
  BluePatternPulseCount=2; BluePatternPulseMs=5; BluePatternPeriodMs=20;
  sendByte('f'); assert(pins[ShutterBlue]==LOW);
  advance(5); assert(pins[ShutterBlue]==HIGH && behaviorBlue.remaining==1);
  advance(15); assert(pins[ShutterBlue]==LOW);
  advance(5); assert(pins[ShutterBlue]==HIGH && !behaviorBlue.active);
  sendByte('F'); assert(behaviorRed.active && pins[ShutterRed]==LOW);
  sendByte('R'); advance(RedPatternPeriodMs); assert(!behaviorRed.active && pins[ShutterRed]==HIGH);
  // Invalid patterns never open the shutter.
  BluePatternPeriodMs=1;
  sendByte('f'); assert(!behaviorBlue.active && pins[ShutterBlue]==HIGH);
  BluePatternPeriodMs=20;
  // End is an immediate stop, preserving every recorded session count.
  std::vector<int*> counts={&TrialNum,&PositiveNum,&NegativeNum,&ManualPositiveNum,&ManualNegativeNum};
#if BEHAVIOR_LEGACY_REWARD_BIG
  counts.push_back(&toneNum);
#else
  counts.insert(counts.end(), {&LeftCueNum,&RightCueNum,&LickCount,&Hit,&Miss,&FalseAlarm,&CorrectReject,&GoNum,&NoGoNum});
#endif
  for (unsigned i=0;i<counts.size();++i) *counts[i]=11+i;
  sendByte('f'); sendByte('F'); sendByte('t');
  // Begin a trial, then emulate in-progress trial outputs.
  sendByte('s');
  digitalWrite(WaterSpout2,HIGH); digitalWrite(Airpuff,HIGH);
  digitalWrite(ShutterBlue,LOW); digitalWrite(ShutterRed,LOW); toneOn=true;
  sendByte('x'); assert(state==Idle); assertClosed();
  advance(10000); sendByte('x'); assertClosed();
  for (unsigned i=0;i<counts.size();++i) assert(*counts[i]==int(11+i));
  // Starting again does not reset recorded counts, either.
  sendByte('s'); assert(state!=Idle);
  Serial.output.clear(); sendByte('?');
  assert(Serial.output.find("BEHAVIOR_STATE task=1 blue=0 red=0")!=std::string::npos);
  for (unsigned i=0;i<counts.size();++i) assert(*counts[i]==int(11+i));
  sendByte('9'); assert(state==Idle); assertClosed();
  // Calibration uses scheduled transitions and can stop while the valve is open.
  sendByte('c'); assert(behaviorCalibrating && behaviorCalibrationDelivered==0);
  advance(CalibrationIntervalMs); assert(pins[WaterSpout2]==HIGH);
  advance(UnitRewardSize); assert(behaviorCalibrationDelivered==1 && pins[WaterSpout2]==LOW);
  advance(CalibrationIntervalMs); assert(pins[WaterSpout2]==HIGH);
  sendByte('x'); assert(!behaviorCalibrating && behaviorCalibrationDelivered==1); assertClosed();
  advance(10000); assertClosed();
  for (unsigned i=0;i<counts.size();++i) assert(*counts[i]==int(11+i));
  CalibrationRepeats=2;
  sendByte('7');
  for (int n=0;n<2;++n) { advance(CalibrationIntervalMs); advance(UnitRewardSize); }
  assert(!behaviorCalibrating && behaviorCalibrationDelivered==2);
  // Manual stimuli cannot collide with a running trial.
  sendByte('s'); sendByte('w');
  assert(!behaviorWater.active && Serial.output.find("ERR BUSY")!=std::string::npos);
  sendByte('x'); assertClosed();
  assert(Serial.output.find("ACK END_TASK counts_preserved")!=std::string::npos);
}
