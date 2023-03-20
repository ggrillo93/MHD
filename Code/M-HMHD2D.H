#include "xtensor.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>

typedef xt::xtensor_fixed<double, xt::xshape<9>> State;
typedef xt::xtensor_fixed<double, xt::xshape<3>> Vec3;
typedef xt::xtensor_fixed<double, xt::xshape<2>> Vec2;
typedef xt::xtensor<double, 3> Tens;
typedef xt::xtensor<double, 2> Mat;
typedef xt::xtensor<double, 1> Arr;

double pressure(State uVec);
double energy(State wVec);
void flux(State uVec, State &fluxVec, unsigned int coord);
void primitive(State uVec, State& wVec);
void conservative(State wVec, State& uVec);
void setBound(Tens& tens);
double fastSpeed(State wVec, unsigned int coord);
void calcWaveSpeeds(State wL, State wR, Vec3& waveSpeeds, unsigned int coord);
void calcUStarK(State wK, State wHLL, double EK, double SK, double vxStar, State& uStarK, unsigned int coord);
void HLLCFlux(State uL, State uR, State& fluxVec, Vec3& waveSpeeds, unsigned int coord);
void computeNetFluxX(Tens uTens, Tens& fluxTens, double step);
void computeNetFluxY(Tens uTens, Tens& fluxTens, double step);
void reconstruct(State uL, State u0, State uR, State& u0Re_L, State& u0Re_R, unsigned int coord);
void halfTimeEvol(State& u0Re_L, State& u0Re_R, double step, unsigned int coord);
double calcMaxSpeed(Tens wTens);
void cleanDiv(State& uL, State& uR, unsigned int coord);
void evolve(Tens& uTens, Tens& fluxTens, double step, unsigned int coord);
void output(std::ofstream& fp, Tens wTens, double t, std::vector<double>& tArr);
void calcDivB(Tens wTens, double dx, std::vector<double>& divBNorm, std::vector<double>& divBMax);
void solveSource(Tens& uTens, double dt);
void updatewTens(Tens uTens, Tens& wTens);
int main();