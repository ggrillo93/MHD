#include "M-HMHD2D.H"
#include "2Dtests.H"
#include "limiters.H"

double Gamma;
unsigned int NGhost = 2, NCellsX, NCellsY, NCCellsX, i, j, k;
double CFL = 0.8, ch, ch2;
const double Cp2ChRatio = 0.18;
std::function<double(double)> limFunc;
State varLim;
std::array<char, 2> BC;
bool divClean;

/* Conservative vector of variables = uVec = (density = rho, rho * vx, rho * vy, rho * vz, total energy = U, Bx, By, Bz, psi)
   Primitive vector of variables = wVec = (rho, vx, vy, vz, pressure = p, Bx, By, Bz, psi) */

double pressure(State uVec) {
    /* Given a conservative state vector uVec, returns the gas pressure */
    double u1 = uVec(0), u2 = uVec(1), u3 = uVec(2), u4 = uVec(3), u5 = uVec(4), u6 = uVec(5), u7 = uVec(6), u8 = uVec(7);
    return (Gamma - 1) * (u5 - 0.5 * ((u2 * u2 + u3 * u3 + u4 * u4) / u1 + u6 * u6 + u7 * u7 + u8 * u8));
}

double energy(State wVec) {
    /* Given a primitive state vector wVec, returns the total energy = U = internal energy + kinetic energy + magnetic energy */
    double rho = wVec(0), vx = wVec(1), vy = wVec(2), vz = wVec(3), p = wVec(4), Bx = wVec(5), By = wVec(6), Bz = wVec(7);
    double v2 = vx * vx + vy * vy + vz * vz;
    double B2 = Bx * Bx + By * By + Bz * Bz;
    return p / (Gamma - 1) + 0.5 * (rho * v2 + B2);
}

void flux(State uVec, State &fluxVec, unsigned int coord) {
    /* Calculates the flux vector along a given direction */
    double u1 = uVec(0), u2 = uVec(1), u3 = uVec(2), u4 = uVec(3), u5 = uVec(4), u6 = uVec(5), u7 = uVec(6), u8 = uVec(7);
    double p = pressure(uVec);
    double B2 = u6 * u6 + u7 * u7 + u8 * u8;
    double pT = p + 0.5 * B2;
    unsigned int vcoord = 1 + coord, BCoord = 5 + coord;
    fluxVec(0) = uVec(vcoord);
    fluxVec(1) = u2 * uVec(vcoord) / u1 - u6 * uVec(BCoord);
    fluxVec(2) = u3 * uVec(vcoord) / u1 - u7 * uVec(BCoord);
    fluxVec(vcoord) += pT;
    fluxVec(3) = u4 * uVec(vcoord) / u1 - u8 * uVec(BCoord);
    fluxVec(4) = (uVec(vcoord) * (u5 + pT) - uVec(BCoord) * (u2 * u6 + u3 * u7 + u4 * u8)) / u1;
    fluxVec(BCoord) = 0.;
    fluxVec(6 - coord) = (u2 * u7 - u6 * u3) / u1 * pow(-1, coord);
    fluxVec(7) = (uVec(vcoord) * u8 - uVec(BCoord) * u4) / u1;
    fluxVec(8) = 0.;
}

void primitive(State uVec, State& wVec) {
    /* Transforms from conservative to primitive variables */
    double u1 = uVec(0);
    double p = pressure(uVec);
    wVec = uVec;
    wVec(1) /= u1;
    wVec(2) /= u1;
    wVec(3) /= u1;
    wVec(4) = p;
}

void conservative(State wVec, State& uVec) {
    /* Transforms from primitive to conservative variables */
    double rho = wVec(0);
    uVec = wVec;
    uVec(1) *= rho;
    uVec(2) *= rho;
    uVec(3) *= rho;
    uVec(4) = energy(wVec);
}

void setBound(Tens& tens) {
    /* Sets the boundary conditions */
    if (BC[0] == 'T' || BC[0] == 'R') {
        State ref = xt::ones<double>({9});
        if (BC[0] == 'R') {
            ref(1) = ref(5) = -1.;
        }
        for (j = NGhost; j < NCellsY - NGhost; j++) { // have not tested with NGhost > 2
            for (i = 0; i < NGhost; i++) {
                xt::view(tens, j, i) = xt::view(tens, j, NGhost * 2 - 1 - i) * ref;
                xt::view(tens, j, NCellsX - NGhost + i) = xt::view(tens, j, NCellsX - NGhost * 2 - i + 1) * ref;
            }
        }
    }
    else if (BC[0] == 'P') {
        for (j = NGhost; j < NCellsY - NGhost; j++) {
            for (i = 0; i < NGhost; i++) {
                xt::view(tens, j, i) = xt::view(tens, j, NCellsX - NGhost * 2 + i);
                xt::view(tens, j, NCellsX - NGhost + i) = xt::view(tens, j, NGhost + i);
            }
        }
    }
    else {
        std::cerr << "BC not implemented";
        std::exit(EXIT_FAILURE);
    }
    if (BC[1] == 'T' || BC[1] == 'R') {
        State ref = xt::ones<double>({9});
        if (BC[1] == 'R') {
            ref(2) = ref(6) = -1.;
        }
        for (i = NGhost; i < NCellsX - NGhost; i++) {
            for (j = 0; j < NGhost; j++) {
                xt::view(tens, j, i) = xt::view(tens, NGhost * 2 - 1 + j, i) * ref;
                xt::view(tens, NCellsY - NGhost + j, i) = xt::view(tens, NCellsY - NGhost * 2 - j + 1, i) * ref;
            }
        }
    }
    else if (BC[1] == 'P') {
        for (i = NGhost; i < NCellsX - NGhost; i++) {
            for (j = 0; j < NGhost; j++) {
                xt::view(tens, j, i) = xt::view(tens, NCellsY - NGhost * 2 + j, i);
                xt::view(tens, NCellsY - NGhost + j, i) = xt::view(tens, NGhost + j, i);
            }
        }
    }
    else {
        std::cerr << "BC not implemented";
        std::exit(EXIT_FAILURE);
    }
}

double fastSpeed(State wVec, unsigned int coord) {
    /* Calculates the MHD fast speed */
    double rho = wVec(0), p = wVec(4), Bx = wVec(5), By = wVec(6), Bz = wVec(7);
    double B2 = Bx * Bx + By * By + Bz * Bz;
    double c_a2 = B2 / rho;
    double c_s2 = Gamma * p / rho;
    double speedSum = c_a2 + c_s2;
    double Bn = wVec(5 + coord);
    return sqrt(0.5 * (speedSum + sqrt(speedSum * speedSum - 4 * c_s2 * Bn * Bn / rho)));
}

void calcWaveSpeeds(State wL, State wR, Vec3& waveSpeeds, unsigned int coord) {
    /* Calculates the wavespeed estimates to be used in the HLLC solver */
    double SL, SR, vStar, maxFast;
    double rhoR = wR(0), vR = wR(1 + coord), pR = wR(4), BxR = wR(5), ByR = wR(6), BzR = wR(7);
    double rhoL = wL(0), vL = wL(1 + coord), pL = wL(4), BxL = wL(5), ByL = wL(6), BzL = wL(7);
    double BnR = wR(5 + coord);
    double BnL = wL(5 + coord);
    double B2R = BxR * BxR + ByR * ByR + BzR * BzR;
    double B2L = BxL * BxL + ByL * ByL + BzL * BzL;
    double pModL = pL + 0.5 * B2L - BnL * BnL;
    double pModR = pR + 0.5 * B2R - BnR * BnR;
    maxFast = std::max(fastSpeed(wL, coord), fastSpeed(wR, coord));
    SL = std::min(vL, vR) - maxFast;
    SR = std::max(vL, vR) + maxFast;
    vStar = (rhoR * vR * (SR - vR) - rhoL * vL * (SL - vL) + pModL - pModR) / (rhoR * (SR - vR) - rhoL * (SL - vL));
    waveSpeeds(0) = SL;
    waveSpeeds(1) = vStar;
    waveSpeeds(2) = SR;
}

void calcUStarK(State wK, State wHLL, double EK, double SK, double vStar, State& uStarK, unsigned int coord) {
    /* Calculates the intermediate state for the HLLC solver */
    double rhoK = wK(0), vxK = wK(1), vyK = wK(2), vzK = wK(3), pK = wK(4), BxK = wK(5), ByK = wK(6), BzK = wK(7);
    double vnK = wK(1 + coord);
    double rhoStarK = rhoK * (SK - vnK) / (SK - vStar);
    double BnK = wK(5 + coord);
    double BnStar;
    if (divClean) {
        BnStar = BnK;
    }
    else {
        BnStar = wHLL(5 + coord);
    }
    double pKT = pK + 0.5 * (BxK * BxK + ByK * ByK + BzK * BzK);
    double pStarT = rhoK * (SK - vnK) * (vStar - vnK) + BnStar * BnStar + pKT - BnK * BnK;
    Vec3 vHLL = {wHLL(1), wHLL(2), wHLL(3)}, BHLL = {wHLL(5), wHLL(6), wHLL(7)};
    double vdotBK = vxK * BxK + vyK * ByK + vzK * BzK;
    double vdotBHLL = xt::sum(BHLL * vHLL)(0);
    uStarK(0) = rhoStarK;
    uStarK(1 + coord) = rhoStarK * vStar;
    uStarK(2 - coord) = rhoStarK * wK(2 - coord) - (BnStar * BHLL(1 - coord) - BxK * ByK) / (SK - vStar);
    uStarK(3) = rhoStarK * vzK - (BnStar * BHLL(2) - BnK * BzK) / (SK - vStar);
    uStarK(4) = EK * (rhoStarK / rhoK) + ((pStarT * vStar - pKT * vnK) - (BnStar * vdotBHLL - BnK * vdotBK)) / (SK - vStar);
    uStarK(5 + coord) = BnStar;
    uStarK(6 - coord) = BHLL(1 - coord);
    uStarK(7) = BHLL(2);
    uStarK(8) = 0.;
}

void HLLCFlux(State uL, State uR, State& fluxVec, Vec3& waveSpeeds, unsigned int coord) {
    /* Solves the Riemann problem between states uL and uR */
    double SL, SR, vStar;
    State uStarK, fluxL, fluxR, wL, wR, uHLL, wHLL;
    primitive(uL, wL);
    primitive(uR, wR);
    calcWaveSpeeds(wL, wR, waveSpeeds, coord);
    SL = waveSpeeds(0), vStar = waveSpeeds(1), SR = waveSpeeds(2);
    if (SL > 0.) {
        flux(uL, fluxVec, coord);
    }
    else if ((vStar >= 0) || (vStar <= 0. && SR >= 0.)) {
        flux(uL, fluxL, coord);
        flux(uR, fluxR, coord);
        uHLL = (SR * uR - SL * uL - fluxR + fluxL) / (SR - SL);
        primitive(uHLL, wHLL);
        if (vStar >= 0) {
            calcUStarK(wL, wHLL, uL(4), SL, vStar, uStarK, coord);
            fluxVec = fluxL + SL * (uStarK - uL);
        }
        else {
            calcUStarK(wR, wHLL, uR(4), SR, vStar, uStarK, coord);
            fluxVec = fluxR + SR * (uStarK - uR);
        }
    }
    else {
        flux(uR, fluxVec, coord);
    }
    if (divClean) {
        fluxVec(5 + coord) = uL(8);
        fluxVec(8) = ch2 * uL(5 + coord);
    }
}

void computeNetFluxX(Tens uTens, Tens& fluxTens, double step) {
    /* Calculates the flux multidimensional matrix along the x-direction */
    State uLL, uL, u0, uR, uLRe_L, uLRe_R, u0Re_L, u0Re_R, fluxR, fluxL;
    Vec3 waveSpeeds;
    for (j = NGhost; j < NCellsY - NGhost; j++) {
        uLL = xt::view(uTens, j, 0);
        uL = xt::view(uTens, j, 1);
        u0 = xt::view(uTens, j, 2);
        uR = xt::view(uTens, j, 3);
        reconstruct(uLL, uL, u0, uLRe_L, uLRe_R, 0);
        halfTimeEvol(uLRe_L, uLRe_R, step, 0);
        reconstruct(uL, u0, uR, u0Re_L, u0Re_R, 0);
        halfTimeEvol(u0Re_L, u0Re_R, step, 0);
        if (divClean) {
            cleanDiv(uLRe_R, u0Re_L, 0);
        }
        HLLCFlux(uLRe_R, u0Re_L, fluxL, waveSpeeds, 0);
        for (i = NGhost; i < NCellsX - NGhost; i++) {
            uLRe_L = u0Re_L;
            uLRe_R = u0Re_R;
            uL = u0;
            u0 = uR;
            uR = xt::view(uTens, j, i + NGhost);
            reconstruct(uL, u0, uR, u0Re_L, u0Re_R, 0);
            halfTimeEvol(u0Re_L, u0Re_R, step, 0);
            if (divClean) {
                cleanDiv(uLRe_R, u0Re_L, 0);
            }
            HLLCFlux(uLRe_R, u0Re_L, fluxR, waveSpeeds, 0);
            xt::view(fluxTens, j, i) = fluxR - fluxL;
            fluxL = fluxR;
        }
    }
}

void computeNetFluxY(Tens uTens, Tens& fluxTens, double step) {
    /* Calculates the flux multidimensional matrix along the y-direction */
    State uDD, uD, u0, uU, uDRe_D, uDRe_U, u0Re_D, u0Re_U, fluxU, fluxD;
    Vec3 waveSpeeds;
    for (i = NGhost; i < NCellsX - NGhost; i++) {
        uDD = xt::view(uTens, NCellsY - 1, i);
        uD = xt::view(uTens, NCellsY - 2, i);
        u0 = xt::view(uTens, NCellsY - 3, i);
        uU = xt::view(uTens, NCellsY - 4, i);
        reconstruct(uDD, uD, u0, uDRe_D, uDRe_U, 1);
        halfTimeEvol(uDRe_D, uDRe_U, step, 1);
        reconstruct(uD, u0, uU, u0Re_D, u0Re_U, 1);
        halfTimeEvol(u0Re_D, u0Re_U, step, 1);
        if (divClean) {
            cleanDiv(uDRe_U, u0Re_D, 1);
        }
        HLLCFlux(uDRe_U, u0Re_D, fluxD, waveSpeeds, 1);
        for (j = NCellsY - NGhost - 1; j > NGhost - 1; j--) {
            uDRe_D = u0Re_D;
            uDRe_U = u0Re_U;
            uD = u0;
            u0 = uU;
            uU = xt::view(uTens, j - NGhost, i);
            reconstruct(uD, u0, uU, u0Re_D, u0Re_U, 1);
            halfTimeEvol(u0Re_D, u0Re_U, step, 1);
            if (divClean) {
                cleanDiv(uDRe_U, u0Re_D, 1);
            }
            HLLCFlux(uDRe_U, u0Re_D, fluxU, waveSpeeds, 1);
            xt::view(fluxTens, j, i) = fluxU - fluxD;
            fluxD = fluxU;
        }
    }
}

void reconstruct(State uL, State u0, State uR, State& u0Re_L, State& u0Re_R, unsigned int coord) {
    /* Linear reconstructs u0 at its left and right boundaries, part of the MUSCL-Hancock scheme */
    double r, limDeltaL, limDeltaR;
    State uDelta, xi;
    for (int ii = 0; ii < u0.size(); ii++) {
        limDeltaR = uR(varLim(ii)) - u0(varLim(ii));
        if (xt::allclose(limDeltaR, 0.)) {
            xi(ii) = 0.;
        }
        else {
            limDeltaL = u0(varLim(ii)) - uL(varLim(ii));
            r = limDeltaL / limDeltaR;
            xi(ii) = limFunc(r);
        }
    }
    uDelta = 0.25 * xi * (uR - uL);
    if (divClean) {
        uDelta(8) = uDelta(5 + coord) = 0.;
    }
    u0Re_L = u0 - uDelta;
    u0Re_R = u0 + uDelta;
}

void halfTimeEvol(State& u0Re_L, State& u0Re_R, double step, unsigned int coord) {
    /* Performs a half time local evolution according to the MUSCL-Hancock scheme */
    State fluxR, fluxL, correction;
    flux(u0Re_L, fluxL, coord);
    flux(u0Re_R, fluxR, coord);
    correction = 0.5 * step * (fluxR - fluxL);
    u0Re_L -= correction;
    u0Re_R -= correction;
}

double calcMaxSpeed(Tens wTens) {
    /* Calculates the maximum wave speed within the computational domain */
    State w0;
    Vec3 maxSpeeds;
    double maxSpeed = 0.;
    for (i = NGhost; i < NCellsX - NGhost; i++) {
        for (j = NGhost; j < NCellsY - NGhost; j++) {
            w0 = xt::view(wTens, j, i);
            maxSpeeds = {fastSpeed(w0, 0) + abs(w0(1)), fastSpeed(w0, 1) + abs(w0(2)), fastSpeed(w0, 2) + abs(w0(3))};
            maxSpeed = std::max(xt::amax(maxSpeeds)(0), maxSpeed);
        }
    }
    return maxSpeed;
}

void cleanDiv(State& uL, State& uR, unsigned int coord) {
    /* Solves the decoupled Riemann problem for Bnorm and psi, used as part of the divergence cleaning scheme */
    unsigned int BCoord = 5 + coord;
    double BnR, BnL, psiR, psiL, BnTilde, psiTilde;
    BnR = uR(BCoord), BnL = uL(BCoord);
    psiR = uR(8), psiL = uL(8);
    BnTilde = 0.5 * (BnL + BnR - (psiR - psiL) / ch);
    psiTilde = 0.5 * (psiR + psiL - ch * (BnR - BnL));
    uR(BCoord) = uL(BCoord) = BnTilde;
    uR(8) = uL(8) = psiTilde;
}

void evolve(Tens& uTens, Tens& fluxTens, double step, unsigned int coord) {
    /* Evolves the multidimensional matrix of states uTens along a given axis */
    if (coord == 0) {
        computeNetFluxX(uTens, fluxTens, step);
    }
    else {
        computeNetFluxY(uTens, fluxTens, step);
    }
    uTens -= step * fluxTens;
}

void output(std::ofstream& fp, Tens wTens, double t, std::vector<double>& tOutArr) {
    /* Outputs quantities to file */
    tOutArr.push_back(t);
    for (k = 0; k < 8; k++) {
        for (j = NGhost; j < NCellsY - NGhost; j++) {
            for (i = NGhost; i < NCellsX - NGhost; i++) {
                fp << wTens(j, i, k) << " ";
            }
            fp << std::endl;
        }
    }
}

void calcDivB(Tens wTens, double dx, std::vector<double>& divBNorm, std::vector<double>& divBMax) {
    /* Calculates the divergence error in the computational domain according to the L1 and Linf norms */
    double divBx, divBy, divB0 = 0, maxDivB = 0;
    double BxR, BxL, ByU, ByD;
    for (j = NGhost; j < NCellsY - NGhost; j++) {
        for (i = NGhost; i < NCellsX - NGhost; i++) {
            BxR = wTens(j, i + 1, 5);
            BxL = wTens(j, i - 1, 5);
            ByU = wTens(j - 1, i, 6);
            ByD = wTens(j + 1, i, 6);
            divBx = 0.5 * (BxR - BxL) / dx;
            divBy = 0.5 * (ByU - ByD) / dx;
            divB0 += abs(divBx + divBy);
            maxDivB = std::max(maxDivB, abs(divBx + divBy));
        }
    }
    divBNorm.push_back(divB0 / (NCCellsX * (NCellsY - NGhost * 2)));
    divBMax.push_back(maxDivB);
}

void solveSource(Tens& uTens, double dt) {
    /* Solves the ODE for the mixed divergence cleaning scheme */
    for (i = NGhost; i < NCellsX - NGhost; i++) {
        for (j = NGhost; j < NCellsY - NGhost; j++) {
            uTens(j, i, 8) *= exp(-dt * ch / Cp2ChRatio);
        }
    }
}

void updatewTens(Tens uTens, Tens& wTens) {
    /* Updates the multidimensional matrix of primitive variables wTens at the end of a timestep */
    State w0;
    for (i = NGhost; i < NCellsX - NGhost; i++) {
        for (j = NGhost; j < NCellsY - NGhost; j++) {
            primitive(xt::view(uTens, j, i), w0);
            xt::view(wTens, j, i) = w0;
        }
    }
}

int main() {
    limFunc = minbee;
    
    std::string tName, divCleanResponse, variableLim;
    unsigned int tNum, nSteps = 0, outSteps;

    std::cout << "Enter test number: ";
    std::cin >> tNum;

    std::cout << "Enter number of x-cells: ";
    std::cin >> NCCellsX;
    NCellsX = NCCellsX + NGhost * 2;

    std::cout << "Enter output frequency: ";
    std::cin >> outSteps;

    std::cout << "Use divergence cleaning? (H = Hyperbolic, M = Mixed, N = None) ";
    std::cin >> divCleanResponse;

    std::cout << "Variable limiting (E = Energy, I = Itself) ";
    std::cin >> variableLim;

    if (divCleanResponse == "H" || divCleanResponse == "M") {
        divClean = true;
    }
    else {
        divClean = false;
    }

    if (variableLim == "E") {
        varLim = xt::ones<int>({9}) * 4;
    }
    else {
        varLim = xt::arange<int>(9);
    }

    double dt, dx, T, step, maxSpeed;

    std::ofstream fp;

    Tens uTens, wTens, fluxTens;

    testSetup(wTens, uTens, T, dx, Gamma, NCellsX, NCellsY, NGhost, BC, tNum, tName);

    fluxTens = xt::zeros<double>(wTens.shape());

    std::vector<double> tOutArr, divBNorm, divBMax, tArr;

    tArr.push_back(0.);
    calcDivB(wTens, dx, divBNorm, divBMax);

    double t = 0;
    
    fp.open("MHDM-H" + tName + divCleanResponse + std::to_string(NCCellsX) + variableLim + ".out");

    output(fp, wTens, t, tOutArr);
    int maxSteps = 10000;

    while (t < T && nSteps < maxSteps) {
        nSteps += 1;
        maxSpeed = calcMaxSpeed(wTens);
        dt = std::min(CFL * dx / maxSpeed, T - t);
        if (divClean) {
            ch = maxSpeed;
            ch2 = ch * ch;
        }
        t += dt;
        std::cout << "Step = " << nSteps << ", dt = " << dt << ", t = " << t << std::endl;
        step = dt/dx;
        evolve(uTens, fluxTens, step, 0); // evolve along the x-direction
        setBound(uTens);
        evolve(uTens, fluxTens, step, 1); // evolve along the y-direction
        if (divCleanResponse == "M") {
            solveSource(uTens, dt);
        }
        updatewTens(uTens, wTens);
        setBound(uTens); // set boundary conditions
        setBound(wTens);
        calcDivB(wTens, dx, divBNorm, divBMax); // calculate divergence error
        tArr.push_back(t);
        if (nSteps % outSteps == 0 || xt::allclose(t, T)) {
            output(fp, wTens, t, tOutArr); // output to file
        }
    }

    for (i = 0; i < tOutArr.size(); i++) {
        fp << tOutArr[i] << " "; // write output times to file
    }

    fp << std::endl;

    for (i = 0; i < tArr.size(); i++) {
        fp << tArr[i] << " "; // write all simulation times to file
    }

    fp << std::endl;

    for (i = 0; i < divBNorm.size(); i++) {
        fp << divBNorm[i] << " "; // write L1 norm of the divergence at each step to file
    }

    fp << std::endl;

    for (i = 0; i < divBMax.size(); i++) {
        fp << divBMax[i] << " "; // write Linf norm of the divergence to file
    }

    fp.close();

    std::cout << "Number of steps = " << nSteps << std::endl;

    return 0;
}