#include <WLMC.h>
#include <armadillo>
#include <vector>
#include <algorithm>
#include <sys/time.h>

#include "ising2d.h"

using namespace arma;
using namespace std;


int main()
{
    srand(time(NULL));

    uint nbins = 1000;

    uint NX = 30;
    uint NY = 30;

    uint Np = NX*NY/2;

    uint overlap = 20;
    uint minWindowSize = 100;

    ising2D system(NX, NY, Np, overlap, minWindowSize);

    WLMC::Window *mainWindow = system.execute(nbins, 0, datum::e, 1 + 1E-6);

    vec E = mainWindow->energies();
    vec DOS = mainWindow->DOS();

    double T, dT;
    uint N = 100000;

    vec temperatures = linspace(0.1, 40.0, N);
    dT = temperatures(1) - temperatures(0);

    vec avgEs(N);

    double partitionFunc, avgE, pdf;
    for (uint i = 0; i < N; ++i)
    {
        T = temperatures(i);

        partitionFunc = 0;
        avgE = 0;

        for (uint bin = 0; bin < nbins; ++bin)
        {
            if (mainWindow->isDeflatedBin(bin))
            {
                continue;
            }

            pdf = DOS(bin)*exp(-E(bin)/T);

            partitionFunc += pdf;
            avgE += pdf*E(bin);

        }

        avgE /= partitionFunc;

        avgEs(i) = avgE;

    }

    vec C(N, fill::zeros);
    mat CT(N, 2);

    CT.col(1) = temperatures;


    for (uint i = 1; i < N - 1; ++i)
    {
        C(i) = (avgEs(i + 1) - avgEs(i - 1))/(2*dT);
    }

    CT.col(0) = C;

    CT.save("CT.arma");

    cout << C.max()/mean(C)
         << " at "
         << temperatures(find(C == C.max()).eval()(0))
         << endl;

    return 0;
}
