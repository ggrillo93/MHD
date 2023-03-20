#include "2Dtests.H"

void setupTensRiemann(Tens& wTens, Tens& uTens, State wR, State wL, Arr xCoord, Arr yCoord, unsigned int NCCellsX, unsigned int NCCellsY, unsigned int NGhost, std::function<bool(double, double)> separator) {
    /* Set up initial multidimensional matrices for Riemann-like initial conditions */
    State uL, uR;
    double x, y;
    conservative(wL, uL);
    conservative(wR, uR);
    for (unsigned int i = 0; i < NCCellsX; i++) {
        x = xCoord(i);
        for (unsigned int j = 0; j < NCCellsY; j++) {
            y = yCoord(j);
            if (separator(x, y)) {
                xt::view(uTens, j+NGhost, i+NGhost) = uR;
                xt::view(wTens, j+NGhost, i+NGhost) = wR;
            }
            else {
                xt::view(uTens, j+NGhost, i+NGhost) = uL;
                xt::view(wTens, j+NGhost, i+NGhost) = wL;
            }
        }
    }
    setBound(uTens);
    setBound(wTens);
}

void setupTensOT(Tens& wTens, Tens& uTens, Arr xCoord, Arr yCoord, unsigned int NCCellsX, unsigned int NCCellsY, unsigned int NGhost, double Gamma) {
    /* Setup initial multidimensional matrices for the Orszag-Tang vortex */
    State wVec, uVec;
    double x, y;
    double pi = xt::numeric_constants<double>::PI;
    wVec(0) = Gamma * Gamma;
    wVec(3) = wVec(7) = wVec(8) = 0.;
    wVec(4) = Gamma;
    for (unsigned int i = 0; i < NCCellsX; i++) {
        x = xCoord(i);
        wVec(2) = sin(2 * pi * x);
        wVec(6) = sin(4 * pi * x);
        for (unsigned int j = 0; j < NCCellsY; j++) {
            y = yCoord(j);
            wVec(1) = -sin(2 * pi * y);
            wVec(5) = wVec(1);
            conservative(wVec, uVec);
            xt::view(uTens, j+NGhost, i+NGhost) = uVec;
            xt::view(wTens, j+NGhost, i+NGhost) = wVec;
        }
    }
    setBound(uTens);
    setBound(wTens);
}

void setupTensKH(Tens& wTens, Tens& uTens, Arr xCoord, Arr yCoord, unsigned int NCCellsX, unsigned int NCCellsY, unsigned int NGhost, double Gamma) {
    /* Setup initial multidimensional matrices for the Kelvin-Helmholtz instability */
    State wVec, uVec;
    double x, y;
    double pi = xt::numeric_constants<double>::PI;
    double theta = pi/3.;
    double preF;
    wVec(0) = 1.;
    wVec(3) = wVec(6) = wVec(8) = 0.;
    wVec(5) = 0.1 * cos(theta);
    wVec(7) = 0.1 * sin(theta);
    wVec(4) = 1./Gamma;
    for (unsigned int i = 0; i < NCCellsX; i++) {
        x = xCoord(i);
        preF = 0.01 * sin(2 * pi * x);
        for (unsigned int j = 0; j < NCCellsY; j++) {
            y = yCoord(j);
            wVec(1) = 0.5 * tanh(20 * y);
            wVec(2) = preF * exp(-y * y / 0.01);
            conservative(wVec, uVec);
            xt::view(uTens, j+NGhost, i+NGhost) = uVec;
            xt::view(wTens, j+NGhost, i+NGhost) = wVec;
        }
    }
    setBound(uTens);
    setBound(wTens);
}

void testSetup(Tens& wTens, Tens& uTens, double& T, double& dx, double& Gamma, unsigned int NCellsX, unsigned int& NCellsY, unsigned int NGhost, std::array<char, 2>& BC, unsigned int tNum, std::string& tName) {
    /* Helper function to setup initial conditions */
    double minX, maxX, minY, maxY;
    unsigned int NCCellsY, NCCellsX = NCellsX - NGhost * 2;
    State wL, wR;
    std::function<bool(double, double)> separator;
    switch(tNum) {
        case(1): { // Toro Test #1
            wL = {1., 0., 0., 0., 1., 0., 0., 0., 0.}, wR = {0.125, 0., 0., 0., 0.1, 0., 0., 0., 0.};
            minX = 0., maxX = 1.;
            minY = minX, maxY = maxX;
            separator = [](double x, double y){return x > 0.5; };
            tName = "Toro1H";
            T = 0.25;
            Gamma = 1.4;
            BC = {'T', 'T'};
            break;
        }
        case(2): {
            wL = {1., 0., 0., 0., 1., 0., 0., 0., 0.}, wR = {0.125, 0., 0., 0., 0.1, 0., 0., 0., 0.};
            minX = 0., maxX = 1.;
            minY = minX, maxY = maxX;
            separator = [](double x, double y){return y > 0.5; };
            tName = "Toro1V";
            T = 0.25;
            Gamma = 1.4;
            BC = {'T', 'T'};
            break;
        }
        case(3): { // Diagonal Toro Test #1
            wL = {1., 0., 0., 0., 1., 0., 0., 0., 0.}, wR = {0.125, 0., 0., 0., 0.1, 0., 0., 0., 0.};
            minX = 0., maxX = 1.;
            minY = minX, maxY = maxX;
            separator = [](double x, double y){return (x + y > 1); };
            tName = "Toro1Diag";
            T = 0.25;
            Gamma = 1.4;
            BC = {'T', 'T'};
            break;
        }
        case(4): { // Cylindrical explosion
            wL = {1., 0., 0., 0., 1., 0., 0., 0., 0.}; // inside
            wR = {0.125, 0., 0., 0., 0.1, 0., 0., 0.}; // outside
            minX = 0., maxX = 2.;
            minY = minX, maxY = maxX;
            separator = [](double x, double y){return (((x - 1) * (x - 1) + (y - 1) * (y - 1)) > 0.16); };
            tName = "CylinderExp";
            T = 0.25;
            Gamma = 1.4;
            BC = {'T', 'T'};
            break;
        }
        case(5): { // Brio and Wu test, horizontal
            wL = {1., 0., 0., 0., 1, 0.75, 1., 0., 0.};
            wR = {0.125, 0., 0., 0., 0.1, 0.75, -1., 0., 0.};
            minX = 0., maxX = 800.;
            dx = (maxX - minX) / NCCellsX;
            minY = minX, maxY = maxX; // can set maxY = 2*dx to make test effectively 1D
            separator = [](double x, double y){return x > 400; };
            tName = "BrioWuX";
            T = 80.;
            Gamma = 2.;
            BC = {'T', 'T'};
            break;
        }
        case(6): { // Brio and Wu test, vertical
            wL = {1., 0., 0., 0., 1, 1., 0.75, 0., 0.};
            wR = {0.125, 0., 0., 0., 0.1, -1., 0.75, 0., 0.};
            minX = 0., maxX = 800.;
            minY = minX, maxY = maxX;
            separator = [](double x, double y){return y > 400; };
            tName = "BrioWuY";
            T = 80.;
            Gamma = 2.;
            BC = {'T', 'T'};
            break;
        }
        case(7): { // Diagonal Brio and Wu test
            double BxL = 0.75, ByL = 1., BxR = 0.75, ByR = -1.;
            wL = {1., 0., 0., 0., 1, (BxL - ByL) / sqrt(2), (BxL + ByL) / sqrt(2), 0., 0.};
            wR = {0.125, 0., 0., 0., 0.1, (BxR - ByR) / sqrt(2), (BxR + ByR) / sqrt(2), 0., 0.};
            minX = 0., maxX = 800.;
            minY = minX, maxY = maxX;
            separator = [](double x, double y){return x + y > 800; };
            tName = "BrioWuDiag";
            T = 80.;
            Gamma = 2.;
            BC = {'T', 'T'};
            break;
        }
        case(8): { // Orszag-Tang, t = 0.5
            Gamma = 5./3.;
            minX = minY = 0.;
            maxX = maxY = 1.;
            tName = "OrszagTangA";
            T = 0.5;
            BC = {'P', 'P'};
            break;
        }
        case(9): { // Orszag-Tang, t = 1
            Gamma = 5./3.;
            minX = minY = 0.;
            maxX = maxY = 1.;
            tName = "OrszagTangB";
            T = 1.;
            BC = {'P', 'P'};
            break;   
        }
        case(10): { // Kelvin-Helmholtz, t = 5
            Gamma = 5./3.;
            minX = 0.;
            minY = -1.;
            maxX = maxY = 1.;
            tName = "KelvinHelmholtzA";
            T = 5.;
            BC = {'P', 'R'};
            break;
        }
        case(11): { // Kelvin-Helmholtz, t = 8
            Gamma = 5./3.;
            minX = 0.;
            minY = -1.;
            maxX = maxY = 1.;
            tName = "KelvinHelmholtzB";
            T = 8.;
            BC = {'P', 'R'};
            break;
        }
        case(12): { // Kelvin-Helmholtz, t = 12
            Gamma = 5./3.;
            minX = 0.;
            minY = -1.;
            maxX = maxY = 1.;
            tName = "KelvinHelmholtzC";
            T = 12.;
            BC = {'P', 'R'};
            break;
        }
        default: {
            std::cerr << "Test not implemented";
            std::exit(EXIT_FAILURE);
            break;
        }
    }
    dx = (maxX - minX) / NCCellsX;
    Arr xCoord = xt::linspace(minX + dx/2., maxX - dx/2., NCCellsX);
    NCCellsY = (maxY - minY) / dx;
    NCellsY = NCCellsY + NGhost * 2;
    std::cout << "y-cells = " << NCCellsY << std::endl;
    Arr yCoord = xt::linspace(maxY - dx/2., minY + dx/2., NCCellsY);
    xt::xarray<int>::shape_type sh = {NCellsY, NCellsX, 9};
    wTens = uTens = xt::zeros<double>(sh);
    if (tNum <= 7) {
        setupTensRiemann(wTens, uTens, wR, wL, xCoord, yCoord, NCCellsX, NCCellsY, NGhost, separator);
    }
    else if (tNum <= 9) {
        setupTensOT(wTens, uTens, xCoord, yCoord, NCCellsX, NCCellsY, NGhost, Gamma);
    }
    else {
        setupTensKH(wTens, uTens, xCoord, yCoord, NCCellsX, NCCellsY, NGhost, Gamma);
    }
}