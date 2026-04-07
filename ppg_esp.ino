#include <Wire.h>
#include "MAX30105.h"

MAX30105 particleSensor;

const long FINGER_THRESHOLD = 50000;

void setup() {
  Serial.begin(115200);
  delay(1000);

  Wire.begin(21, 22);   // SDA = 21, SCL = 22

  if (!particleSensor.begin(Wire, I2C_SPEED_STANDARD)) {
    Serial.println("ERROR");
    while (1);
  }

  // setup(ledBrightness, sampleAverage, ledMode, sampleRate, pulseWidth, adcRange)
  // sampleAverage = 1 keeps effective output rate high
  particleSensor.setup(40, 1, 2, 100, 411, 4096);

  particleSensor.setPulseAmplitudeRed(0x1F);
  particleSensor.setPulseAmplitudeIR(0x1F);
  particleSensor.setPulseAmplitudeGreen(0);

  delay(1000);
}

void loop() {
  particleSensor.check();   // update FIFO

  while (particleSensor.available()) {
    long irValue = particleSensor.getIR();
    unsigned long t = millis();

    Serial.print(t);
    Serial.print(",");

    if (irValue < FINGER_THRESHOLD) {
      Serial.println(0);
    } else {
      Serial.println(irValue);
    }

    particleSensor.nextSample();
  }
}