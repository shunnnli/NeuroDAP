// Exercise real firmware against mocked pins/time, including archived sketches.
void sendByte(char command) { Serial.input.push_back(command); loop(); }
void advance(unsigned long milliseconds) { fakeTime += milliseconds; loop(); }
void assertClosed() {
#if BEHAVIOR_HAS_WATER
  assert(pins[WaterSpout]==LOW && pins[WaterSpout2]==LOW);
  assert(pins[WaterSpout_copy]==LOW && pins[WaterSpout2_copy]==LOW);
#endif
#if BEHAVIOR_HAS_AIR
  assert(pins[Airpuff]==LOW && pins[Airpuff_copy]==LOW);
#endif
#if BEHAVIOR_HAS_TONE
  assert(!toneOn);
#endif
#if BEHAVIOR_HAS_BLUE
  assert(pins[ShutterBlue]==BEHAVIOR_SHUTTER_CLOSED && !behaviorBlue.active);
#endif
#if BEHAVIOR_HAS_RED
  assert(pins[ShutterRed]==BEHAVIOR_SHUTTER_CLOSED && !behaviorRed.active);
#endif
}
int main() {
  // Keep physical calibration buttons and Bonsai's trigger inactive.
  inputOverrides[4]=inputOverrides[11]=inputOverrides[12]=inputOverrides[53]=LOW;
  setup();
  assert(Serial.output.find("BEHAVIOR_CONTROLS 1")!=std::string::npos);
  assert(Serial.output.find("BEHAVIOR_STATE task=0 blue=0 red=0 calibration=0")!=std::string::npos);
  sendByte('x'); assert(state==Idle); assertClosed();
#if BEHAVIOR_HAS_WATER
  SmallRewardSize=11; BigRewardSize=37;
  sendByte('w'); assert(behaviorWater.active && behaviorWater.duration==11 && pins[WaterSpout2]==HIGH);
  advance(11); assert(!behaviorWater.active && pins[WaterSpout2]==LOW);
  sendByte('W'); assert(behaviorWater.duration==37 && pins[WaterSpout2_copy]==HIGH);
  advance(37); assert(!behaviorWater.active && pins[WaterSpout2_copy]==LOW);
  CalibrationRepeats=2;
  sendByte('c'); assert(behaviorCalibrating);
  Serial.output.clear(); sendByte('?');
  assert(Serial.output.find("calibration=1")!=std::string::npos);
  advance(CalibrationIntervalMs); assert(pins[WaterSpout2]==HIGH);
  sendByte('x'); assert(!behaviorCalibrating); assertClosed();
  advance(10000); assertClosed();
  sendByte('c');
  for(int i=0;i<2;++i) { advance(CalibrationIntervalMs); advance(UnitRewardSize); }
  assert(!behaviorCalibrating && behaviorCalibrationDelivered==2);
  assert(Serial.output.find("DONE CALIBRATION deliveries=2")!=std::string::npos);
#else
  for(char c : {'w','W','c'}) {
    Serial.output.clear(); sendByte(c);
    assert(Serial.output.find("ERR UNSUPPORTED")!=std::string::npos);
  }
#endif
#if BEHAVIOR_HAS_AIR
  sendByte('p'); assert(behaviorAir.active && pins[Airpuff]==HIGH);
  advance(SmallPunishSize); assert(!behaviorAir.active && pins[Airpuff]==LOW);
#else
  Serial.output.clear(); sendByte('p'); assert(Serial.output.find("ERR UNSUPPORTED")!=std::string::npos);
#endif
#if BEHAVIOR_HAS_TONE
  sendByte('t'); assert(toneOn);
  advance(ShortToneDuration); assert(!toneOn);
#endif
#if BEHAVIOR_HAS_BLUE
  sendByte('b'); assert(pins[ShutterBlue]==BEHAVIOR_SHUTTER_OPEN);
  Serial.output.clear(); sendByte('?');
  assert(Serial.output.find("blue=1")!=std::string::npos);
  sendByte('B'); assertClosed();
  sendByte('f'); assert(behaviorBlue.active && pins[ShutterBlue]==BEHAVIOR_SHUTTER_OPEN);
  advance(BluePatternPulseMs); assert(pins[ShutterBlue]==BEHAVIOR_SHUTTER_CLOSED);
  sendByte('B'); advance(10000); assertClosed();
#endif
#if BEHAVIOR_HAS_RED
  sendByte('r'); assert(pins[ShutterRed]==BEHAVIOR_SHUTTER_OPEN);
  sendByte('R'); assertClosed();
  sendByte('F'); assert(behaviorRed.active && pins[ShutterRed]==BEHAVIOR_SHUTTER_OPEN);
  advance(RedPatternPulseMs); assert(pins[ShutterRed]==BEHAVIOR_SHUTTER_CLOSED);
  sendByte('R'); advance(10000); assertClosed();
  Serial.output.clear(); sendByte('F');
  for(unsigned long i=0;i<RedPatternPulseCount;++i) {
    advance(RedPatternPulseMs);
    if(i + 1 < RedPatternPulseCount) advance(RedPatternPeriodMs - RedPatternPulseMs);
  }
  assert(!behaviorRed.active && pins[ShutterRed]==BEHAVIOR_SHUTTER_CLOSED);
  assert(Serial.output.find("DONE RED_PATTERN")!=std::string::npos);
#endif
  std::vector<int*> counts=BEHAVIOR_TEST_COUNTS;
  for(unsigned i=0;i<counts.size();++i) *counts[i]=30+i;
#if BEHAVIOR_HAS_TASK
  sendByte('s'); assert(state!=Idle);
#if BEHAVIOR_HAS_WATER
  sendByte('w'); assert(!behaviorWater.active);
#endif
  sendByte('x'); assertClosed();
  advance(10000); assertClosed();
  sendByte('s'); assert(state!=Idle);
  sendByte('9'); assertClosed();
  for(unsigned i=0;i<counts.size();++i) assert(*counts[i]==int(30+i));
#else
  Serial.output.clear(); sendByte('s'); assert(state==Idle);
  assert(Serial.output.find("ERR UNSUPPORTED")!=std::string::npos);
#endif
#ifdef CHECK_LEGACY_BEHAVIOR
  sendByte('1'); assert(behaviorWater.duration==(BEHAVIOR_LEGACY_REWARD_BIG ? BigRewardSize : SmallRewardSize));
  sendByte('x'); sendByte('2'); assert(behaviorAir.active); sendByte('x');
  sendByte('s');
  Serial.input.push_back('x');
  assert(!behaviorWait(500)); assertClosed(); assert(state==Idle);
#endif
#ifdef CHECK_RAMP
  sendByte('s'); Serial.input.push_back('x');
  playRampTone(500,3000,12000); assert(state==Idle); assertClosed();
  sendByte('s'); Serial.input.push_back('x');
  playFluctTone(); assert(state==Idle); assertClosed();
#endif
#ifdef CHECK_BONSAI
  sendByte('s'); inputOverrides[triggerPin]=HIGH;
  advance(100); assert(isStimulating && pins[stimOutPin]==LOW);
  sendByte('x'); assert(!isStimulating); assertClosed();
  advance(1000); assertClosed(); // High trigger must not restart after End task.
#endif
#ifdef CHECK_CALIBRATION_BOX
  sendByte('1'); assert(behaviorCalibrating && behaviorCalibrationValve==WaterSpout);
  advance(100); assert(pins[WaterSpout]==HIGH);
  sendByte('x'); assertClosed(); advance(10000); assertClosed();
  sendByte('2'); assert(behaviorCalibrating && behaviorCalibrationValve==WaterSpout2);
  sendByte('x'); assertClosed();
  sendByte('5'); assert(checkToneSteps && toneOn);
  sendByte('x'); advance(1000); assertClosed(); assert(!checkToneSteps);
#endif
}
