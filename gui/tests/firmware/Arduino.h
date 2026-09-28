// Host-side Arduino mock for execution tests; never drives physical hardware.
#pragma once
#include <cassert>
#include <sstream>
#include <string>
#include <vector>
#include <deque>
using byte = unsigned char;
using boolean = bool;
const int LOW=0, HIGH=1, OUTPUT=1, INPUT=0;
unsigned long fakeTime=100;
int pins[64]={};
bool toneOn=false;
unsigned long millis() { return fakeTime; }
void delay(unsigned long value) { fakeTime += value; }
void pinMode(byte, int) {}
void digitalWrite(byte pin, int value) { pins[pin]=value; }
int digitalRead(byte pin) { return (pin==22 || pin==24) ? pins[pin] : HIGH; }
int analogRead(byte) { return 1; }
void randomSeed(unsigned long) {}
long random(long value) { return value > 0 ? value/2 : 0; }
long random(long start, long end) { return start+(end-start)/2; }
void tone(byte, unsigned int) { toneOn=true; }
void noTone(byte) { toneOn=false; }
struct SerialMock {
  std::deque<char> input;
  std::string output;
  void begin(long) {}
  int available() { return input.size(); }
  int read() { char value=input.front(); input.pop_front(); return value; }
  template<class T> void print(T value) { std::ostringstream text; text<<value; output+=text.str(); }
  template<class T> void println(T value) { print(value); output+='\n'; }
  void println() { output+='\n'; }
} Serial;
