#ifndef WAVEPARAMETER_H
#define WAVEPARAMETER_H
#include "type.h"
#include "flowfield_host_device.h"
#include "model.h"
#include "model_host_device.h"
//*****wave parameters*****//
const Real Ma{0.02};
const Real MeshSize{CS * TimeStep()};
const Real WaveAmplitude{0.02};
const Real WaveLength{1.21};
const Real WavePeriod{1.0};
const Real L0 = WaveAmplitude * 2;
const Real WaterDepth{0.2};
Real SoundofWater{WaveLength / WavePeriod / Ma};
// wave absorption
const Real LinearConst {200 * L0 / SoundofWater};  // unit 1/T
const int ns{10};  // empirical coefficient
Real DampWidth{2.0 * WaveLength};
const int NumberofMesh {(int)(DampWidth / L0 / MeshSize + 0.5)};
//*****wave parameters*****//
#endif