#include "M-HMHD2D.H"

void setupTensRiemann(Tens& wTens, Tens& uTens, State wR, State wL, Arr xCoord, Arr yCoord, unsigned int NCCellsX, unsigned int NCCellsY, unsigned int NGhost, std::function<bool(double, double)> separator);
void setupTensOT(Tens& wTens, Tens& uTens, Arr xCoord, Arr yCoord, unsigned int NCCellsX, unsigned int NCCellsY, unsigned int NGhost, double Gamma);
void setupTensKH(Tens& wTens, Tens& uTens, Arr xCoord, Arr yCoord, unsigned int NCCellsX, unsigned int NCCellsY, unsigned int NGhost, double Gamma);
void testSetup(Tens& wTens, Tens& uTens, double& T, double& dx, double& Gamma, unsigned int NCCellsX, unsigned int& NCCellsY, unsigned int NGhost, std::array<char, 2>& BC, unsigned int tNum, std::string& tName);