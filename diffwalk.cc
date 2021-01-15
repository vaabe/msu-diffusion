#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <time.h>
#include <algorithm>
#include <complex>
#include <random>
#include <cstring>
#include <omp.h>

using namespace std;

double ran1(long *idum);
double diffclock(clock_t clock1, clock_t clock2)
{
    double diffticks = clock1 - clock2;
    double diffms = (diffticks * 10) / CLOCKS_PER_SEC;
    return diffms;
}

// simulation constants
#define D1 1 // diffusion coefficient 1 (in um^2/ms)
#define D2 0.5 // diffusion coefficient 2
#define dt 0.008 // time step in ms
#define abar 6.0 // abar in um

const double P = .4; // permeability in um/ms
const double zeta = 1/abar/P;
const double tau_R = abar * 0.5 / P;
const double TTu = 200 * tau_R; //250 * tau_R;
const int TN = TTu / dt; // number of time steps
const double TT = TN * dt; // total simulation time
const int sub = 50;
const long NP = 1000; // number of random walkers
const int N_m = 1000; // number of membranes
const int rep = 1; // replications with same barriers and initial positions
long iseed, *idum;
const double step1 = sqrt(2 * dt * D1);
const double step2 = sqrt(2 * dt * D2);
const double v1 = step1 / dt;
const double v2 = step2 / dt;
// const double pEX = step * P / D;
double ran;
int regen;

// mersenne twister generator
std::mt19937 generator(iseed);
std::uniform_real_distribution<double> dis(0.0, 1.0); // for random values 0-1
std::uniform_real_distribution<double> epsilon(-abar/2, abar/2); // for random epsilon values (strong disorder)

// generate barriers
void barriergen(double M[], char* barriertype) {
    
    double a;
    double mu = 1.75;
    double amin = (abar*(mu-1))/mu;
    int i = 0, j = 0;

    // periodic
    if (strncmp(barriertype , "periodic", 8) == 0) {
        while (i < N_m - 1) {
            M[i] = j * abar;
            i++;
            j++;
        }
    }
    
    // hyperuniform disorder
    if (strncmp(barriertype, "hyperuniform", 12) == 0) {
        while (i < N_m - 1) {
            M[i] = j*abar + epsilon(generator);
            i++;
            j++;
        }
    }
    
    // strong disorder
    if (strncmp(barriertype, "strong", 6) == 0) {
        while (i < N_m - 1) {
            a += amin * pow(dis(generator), -1/mu);
            M[i] = a;
            i++;
        }
    }
    
    ofstream barriers("./data/barriers.txt");
    
    for (i = 0; i < N_m; i++) {
        barriers << M[i] << endl;
    }
}

int main (int argc, char *argv[]) {
    
    if (strncmp(argv[2], "periodic", 8) != 0 && 
		strncmp(argv[2], "hyperuniform", 12) != 0 && 
		strncmp(argv[2], "strong", 6) != 0) {

        cout << "Argument 2 must be periodic, hyperuniform, or strong" << endl;
        return 999;
    }
    
    if (argv[3] != NULL) {
        regen = atoi(argv[3]);
    }
    else {
        regen = 1;
    }

    // indices
    int i, j, ii, jk, jl, Ntraj;
    clock_t begin = clock();
    clock_t end = clock();
    clock_t begin1 = clock();
    clock_t end1 = clock();
    clock_t begin2 = clock();
    clock_t end2 = clock();
    cout << "Start of simulation... " << endl;
    cout << "zeta " << zeta << endl;
    cout << "Permeability " << P << endl;
    cout << "Number of random walkers " << NP << endl;
    cout << "Time points " << TN << endl;
    cout << "Regenerations " << regen << endl;
    // cout << "pEX " << pEX << endl;
    cout << "step1 " << step1 << endl;
    cout << "step2 " << step2 << endl;
    cout << "v1 " << v1 << endl;
    cout << "v2 " << v2 << endl;
    cout << "tau_R " << tau_R << endl;
    idum = new long;

    iseed = -atoi(argv[1]);
    *idum = -iseed;
    cout << "random seed " << -iseed << endl;

    // simulation variables
    int Tsub = (int)ceil(TN / sub); // sub-sampling
    long double *X, *M22, *M02, *M20, *M13, *M11, *M31, *M40, *Z22, *Z13;
    X = new long double[TN];
    M20 = new long double[Tsub];
    M40 = new long double[Tsub];
    M11 = new long double[(int)Tsub / 2 + 1]; M22 = new long double[(int)Tsub / 2 + 1]; M02 = new long double[(int)Tsub / 2 + 1]; M13 = new long double[(int)Tsub / 2 + 1]; M31 = new long double[(int)Tsub / 2 + 1];
    Z22 = new long double[(int)Tsub / 2 + 1]; Z13 = new long double[(int)Tsub / 2 + 1];
    
    for (ii = 0; ii < TN; ii++) {
        X[ii] = 0;
    }
    for (ii = 0; ii < Tsub; ii++) {
        M20[ii] = 0; M40[ii] = 0;
    }
    for (ii = 0; ii < Tsub / 2; ii++) {
        M11[ii] = 0; M22[ii] = 0; M02[ii] = 0; M13[ii] = 0; M31[ii] = 0;
        Z22[ii] = 0; Z13[ii] = 0;
    }
    
    double xpp, rest, xnp, tru, pex, dt1, dt2, X2;
    long xR;
    int jump, NT, dir, k, ix;
    
    double *xip;
    xip = new double [NP];
    
    long *x_m;
    x_m = new long [NP];
    
    // simulation
    cout << "Starting simulations..." << endl;
    begin = clock();

    for (jl = 0; jl < regen; jl++) {
        
        cout << "Simulation " << jl+1 << ":" << endl;
        begin1 = clock();

        // initialize positions of barriers s
        cout << "Generating barriers... ";
        begin2 = clock();
        
        double M[N_m];
        barriergen(M, argv[2]); // call function to generate barrier positions
        const double Xmax = M[N_m - 1];
        
        end2 = clock();
        cout << "Done! Elapsed time " << double(diffclock(end2, begin2)) << " ms" << endl;

        // initialize positions of random walkers
        cout << "Initializing random points... ";
        begin2 = clock();
        
        for (i = 0; i < NP; i++) { // initial positions of random walkers
            xip[i] = M[0]+dis(generator)*(Xmax-M[0]);
        }
        sort(xip, xip + NP);
        
        // define nearest RH barrier for all points
        j = 0; i = 0;
        while (j < N_m - 1 && i <= (NP - 1)) {
            while (xip[i] < M[j] && i <= NP - 1) {
                x_m[i] = j;
                i++;
            }
            j++;
        }
        for (ii = i; ii < NP; ii++) {
            x_m[ii] = j;
        }
        
        // define step size of particle in each region (alternates between barriers)
        double stepsize[N_m];
        
        for (i = 0; i < N_m - 1; i = i + 2) {
            stepsize[i] = step1;
        }
        for (j = 1; j < N_m - 1; j = j + 2) {
            stepsize[j] = step2;
        }
        
        // define speed of particle in each region (alternates between barriers)
        double v[N_m];
        
        for (i = 0; i < N_m - 1; i = i + 2) {
            v[i] = v1;
        }
        for (j = 1; j < N_m - 1; j = j + 2) {
            v[j] = v2;
        }

        end2 = clock();
        cout << "Done! Elapsed time " << double(diffclock(end2, begin2)) << " ms" << endl;
        
        // start simulating pseudo random walk
        cout << "Simulating pseudo random walk... " << endl;
        begin2 = clock();
        
        for (jk = 0; jk < rep; jk++) {
            for (k = 0; k < NP; k++) { // loop over all random walkers (NP = number of particles)
                // initializing
                xpp = xip[k]; // random positions from the initialization process
                jump = 0;
                NT = 0; // dummy index
                xR = x_m[k]; // nearest RH barrier
                for (tru = dt; tru <= TT; tru = tru + dt) { // simulation for 1 random walker (TT = total time)
                    // cout << "position " << xpp << endl;
                    // cout << "XMAX " << Xmax << endl;
                    dir = (floor(dis(generator) + 0.5)) * 2 - 1; // set direction right/left
                    rest = stepsize[xR] - max(xpp - M[xR + (dir - 1) / 2], M[xR + (dir - 1) / 2] - xpp);
                    if (rest <= 0) { // if don't cross a barrier, move at constant speed
                        xnp = xpp + dir * v[xR] * dt;
                    }
                    else { // if cross a barrier change speed to reflect new diffusion constant
                        dt1 = rest / v[xR];
                        dt2 = dt - dt1;
                        xnp = xpp + dir * (v[xR] * dt1 + v[xR + 1] * dt2);
                        
                        // check for edges (periodic boundary conditions)
                        if (xR == 0) {
                            xR = N_m - 1;
                            jump -= 1;
                            xnp = xnp + Xmax;
                        }
                        else if (xR == N_m) {
                            xR = 1;
                            jump++;
                            xnp = xnp - Xmax;
                        }
                        
                        rest = rest - (M[xR] + M[xR - 1]);
                    }
                    xpp = xnp; // set old position to new position
                    
                    X[NT] = (xpp + jump * Xmax - xip[k]); // displacement
                    if (isnan(xpp)) {
                        cout << "NAN in " << ": NT = " << NT << " " << xnp << " " << jump << " " << " " << dir << " " << Xmax << " " << xip[k] << endl;
                        exit(1);
                    }
                    NT++;
                } // end simulation for 1 random walker (tru)
                
                for (ii = 0; ii < Tsub; ii++) { // update statistics for one particle for each time step (add one for each loop)
                    // M10[i] += X[i];
                    i = ii * sub;
                    M20[ii] += X[i] * X[i];
                    M40[ii] += X[i] * X[i] * X[i] * X[i];
                    
                    if (i < NT / 2) {
                        ix = 2 * i;
                        X2 = X[ix] - X[i];
                        // M01[i]+=X2;
                        M02[ii] += X2 * X2;
                        M11[ii] += X2 * X[i];
                        M22[ii] += X2 * X2*X[i] * X[i];
                        M31[ii] += X[i] * X[i] * X[i] * X2;
                        M13[ii] += X[i] * X2*X2*X2;
                    }
                }
            } // end loop over all particles (k)
        } // end loop over all replications
        end2 = clock();
        cout << "Done! Elapsed time " << double(diffclock(end2, begin2)) / 10 << " s" << endl;
        end1 = clock();
        cout << "End of simulation " << jl+1 << ". Elapsed time " << double(diffclock(end1, begin1)) / 10 << " s" << endl;
    }
    
    // end loop over all regenerations
    
    delete xip; delete x_m;  delete[] X;
    
    Ntraj = NP * rep * regen;
    for (i = 0; i < Tsub; i++) { // calculate the diffusivity/SM/FM
        M20[i] = M20[i] / Ntraj;
        M40[i] = M40[i] / Ntraj;
        if (i < Tsub / 2) {
            Z22[i] = M22[i] / Ntraj - M20[i] * M02[i] / Ntraj - 2 * M11[i] / Ntraj*M11[i] / Ntraj;
            Z13[i] = (M31[i] / Ntraj + M13[i] / Ntraj - 3 * M11[i] / Ntraj*(M02[i] / Ntraj + M20[i])) / 2;
        }
    }
    
    delete[] M02; delete[] M11; delete[] M22;
    delete[] M31; delete[] M13;
    
    end = clock();
    cout << "End of all simulations. Elpased time " << double(diffclock(end, begin)) / 10 << " s" << endl;
    cout << "Writing data to file..." << endl;;
    
    char M2filename[200];
    char M4filename[200];
    char Z22filename[200];
    char timefilename[200];
    char Z13filename[200];
    
    sprintf(M2filename, "./data/M2_%.3g_TT_%.3g_TN_%.3g_NP_%.3g_regen_%d_rep_%d_P_%.3g_zeta%.3g_tauR_%.3g_abar_%.3g", float(N_m), TT, float(TN), float(NP), regen, rep, P, zeta, tau_R, abar);
    sprintf(M4filename, "./data/M4_%.3g_TT_%.3g_TN_%.3g_NP_%.3g_regen_%d_rep_%d_P_%.3g_zeta%.3g_tauR_%.3g_abar_%.3g", float(N_m), TT, float(TN), float(NP), regen, rep, P, zeta, tau_R, abar);
    sprintf(Z22filename, "./data/Z22_%.3g_TT_%.3g_TN_%.3g_NP_%.3g_regen_%d_rep_%d_P_%.3g_zeta%.3g_tauR_%.3g_abar_%.3g", float(N_m), TT, float(TN), float(NP), regen, rep, P, zeta, tau_R, abar);
    sprintf(timefilename, "./data/time_%.3g_TT_%.3g_TN_%.3g_NP_%.3g_regen_%d_rep_%d_P_%.3g_zeta%.3g_tauR_%.3g_abar_%.3g", float(N_m), TT, float(TN), float(NP), regen, rep, P, zeta, tau_R, abar);
    sprintf(Z13filename, "./data/Z13_%.3g_TT_%.3g_TN_%.3g_NP_%.3g_regen_%d_rep_%d_P_%.3g_zeta%.3g_tauR_%.3g_abar_%.3g", float(N_m), TT, float(TN), float(NP), regen, rep, P, zeta, tau_R, abar);
    
    ofstream M2file(M2filename);
    ofstream M4file(M4filename);
    ofstream Z22file(Z22filename);
    ofstream timefile(timefilename);
    ofstream Z13file(Z13filename);
    
    cout << "size = " << TN << std::endl;
    
    M2file << "% NP = " << NP << "\n%regen = " << regen <<  "\n%rep = " << rep << "\n%N_m = " << float(N_m) << "\n%Time steps = " << TN << "\n%dt = " << dt << "\n%P = " << P << "\n%abar = " << abar << "\n%zeta = " << zeta  << endl;
    timefile << "% NP = " << NP << "regen = " << regen << "rep = " << rep << "\n%N_m = " << float(N_m) << "\n%Time steps = " << TN << "\n%dt = " << dt << "\n%P = " << P << "\n%abar = " << abar << "\n%zeta = " << zeta << endl;
    M4file << "% NP = " << NP << "\n%regen = " << regen << "rep = " << rep << "\n%N_m = " << float(N_m) << "\n%Time steps = " << TN << "\n%dt = " << dt << "\n%P = " << P << "\n%abar = " << abar << "\n%zeta = " << zeta << endl;
    Z22file << "% NP = " << NP << "\n%regen = " << regen << "rep = " << rep << "\n%N_m = " << float(N_m) << "\n%Time steps = " << TN << "\n%dt = " << dt << "\n%P = " << P << "\n%abar = " << abar << "\n%zeta = " << zeta << endl;
    Z13file << "% NP = " << NP << "\n%regen = " << regen << "rep = " << rep << "\n%N_m = " << float(N_m) << "\n%Time steps = " << TN << "\n%dt = " << dt << "\n%P = " << P << "\n%abar = " << abar << "\n%zeta = " << zeta << endl;
    
    for (i = 0; i < Tsub; i++) {
        M2file << std::fixed << std::setprecision(8) << M20[i] << endl;
        M4file << std::fixed << std::setprecision(8) << M40[i] << endl;
        timefile << (i)*sub*dt+dt<< endl;
    }
    for (i = 0; i < Tsub / 2; i++) {
        Z22file << std::fixed << std::setprecision(8) << Z22[i] << endl;
        Z13file << std::fixed << std::setprecision(8) << Z13[i] << endl;
    }
    
    M2file.close();
    M4file.close();
    Z22file.close();
    timefile.close();
    Z13file.close();
    delete M20; delete M40;
    delete Z22; delete Z13;
    delete idum;
    cout << "Done!" << endl;
    return 0;
}


#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-12
#define RNMX (1.0-EPS)

double ran1(long *idum) {
    long j;
    long k;
    static long iy = 0;
    static long iv[NTAB];
    double temp;
    
    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum = 1;
        else *idum = -(*idum);
        for (j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ;
            *idum = IA * (*idum - k * IQ) - IR * k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;
    if (*idum < 0) *idum += IM;
    j = iy / NDIV;
    iy = iv[j];
    iv[j] = *idum;
    if ((temp = AM * iy) > RNMX) return RNMX;
    else return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

