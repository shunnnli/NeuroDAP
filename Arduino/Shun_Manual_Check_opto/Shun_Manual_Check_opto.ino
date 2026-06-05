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

  Serial.println("type 1 for left reward / 2 for right reward");
  Serial.println("for laser shutter: type 3 for blue / 4 for red");
}

void loop()
{
  // Start button check
  if (digitalRead(Button) == 1) {
    Serial.println("Start Button pushed");
  }

  // Left SideButtonsForRewards callibration (drop X 50 times)
  if (digitalRead(Button1) == 1 && ButtonStatus1 == 0)
  {
    ButtonStatus1 = 1;

    Serial.println("Left reward button pushed");

    for (int i = 0; i < num_repeat; i++) {
      digitalWrite(WaterSpout, HIGH);
      digitalWrite(WaterSpout_copy, HIGH);
      delay(RewardSizeLeft);
      digitalWrite(WaterSpout, LOW);
      digitalWrite(WaterSpout_copy, LOW);
      delay(50);
    }
    Serial.print("Finished delivery x");
    Serial.println(num_repeat);
    ButtonStatus1 = 0;
  }

  // RIGHT SideButtonsForRewards callibration (drop X 50 times)
  if (digitalRead(Button2) == 1 && ButtonStatus2 == 0)
  {
    ButtonStatus2 = 1;

    Serial.println("Right reward button pushed");

    for (int i = 0; i < num_repeat; i++) {
      digitalWrite(WaterSpout2, HIGH);
      digitalWrite(WaterSpout2_copy, HIGH);
      delay(RewardSizeRight);
      digitalWrite(WaterSpout2, LOW);
      digitalWrite(WaterSpout2_copy, LOW);
      delay(100);
    }
    Serial.print("Finished delivery x");
    Serial.println(num_repeat);
    ButtonStatus2 = 0;
  }

  // Lick Detection //
  // LEFT
  if (digitalRead(LickDetect1) == 0)
  {
    Serial.println("Left lick detected");
  }

  // RIGHT
  if (digitalRead(LickDetect2) == 0)
  {
    Serial.println("Right lick detected");
  }

  // Serial input from user for dispensing //
  if (Serial.available() > 0) {
    // read the incoming byte:
    SerialInput = Serial.read();

    Serial.print("I received: ");
    Serial.println(SerialInput);

    if (SerialInput == '1') // dispense left reward
    {
      for (int i = 0; i < num_repeat; i++) {
        digitalWrite(WaterSpout, HIGH);
        digitalWrite(WaterSpout_copy, HIGH);
        delay(RewardSizeLeft);
        digitalWrite(WaterSpout, LOW);
        digitalWrite(WaterSpout_copy, LOW);
        delay(100);
      }
      Serial.print("Finished delivery x");
      Serial.println(num_repeat);
      //Serial.println(" (unit size:)");
    }

    if (SerialInput == '2') // disepnse right reward
    {
      for (int i = 0; i < num_repeat; i++) {
        digitalWrite(WaterSpout2, HIGH);
        digitalWrite(WaterSpout2_copy, HIGH);
        delay(RewardSizeRight);
        digitalWrite(WaterSpout2, LOW);
        digitalWrite(WaterSpout2_copy, LOW);
        delay(100);
      }
      Serial.print("Finished delivery x");
      Serial.println(num_repeat);
    }

    if (SerialInput == '3') {
      digitalWrite(ShutterBlue, HIGH);
      //      delay(2000);
      //      digitalWrite(ShutterBlue, LOW);
    }

    if (SerialInput == '4') {
      digitalWrite(ShutterRed, HIGH);
      //      digitalWrite(ShutterRed, HIGH);
      //      delay(2000);
      //      digitalWrite(ShutterRed, LOW);
      //
      //      for (int i = 0; i < 100; i++) {
      //        digitalWrite(ShutterRed, HIGH);
      //        delay(2);
      //        digitalWrite(ShutterRed, LOW);
      //        delay(48);
      //      }
    }

    if (SerialInput == '5') {
      // speaker check x 3 times
      for (int i = 0; i < 50; i++) {
        tone(Speaker, GoCueFreq);
        delay(100);
        tone(Speaker, NoGoCueFreq);
        delay(100);
      }
      noTone(Speaker);
    }

    delay(1000); //delay so that reward button is no longer pressed for one drop
  }
}
