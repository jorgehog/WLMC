#include <WLMC.h>
#include <armadillo>
#include <vector>
#include <algorithm>
#include <sys/time.h>

#include "harmosc.h"
//#include "ising2d.h"

#include <DCViz/include/DCViz.h>

#include <unittest++/UnitTest++.h>

#include <iomanip>

using namespace arma;
using namespace std;

const WLMC::Window &harmosc();
//void ising();

TEST(HARMOSC)
{
    const WLMC::Window & mainWindow = harmosc();

    vec compareDOS = mainWindow.logDOS()/mainWindow.logDOS(1);

    double overlap = 0;
    for (uint bin = 2; bin < mainWindow.nbins(); ++bin)
    {
        overlap += fabs(compareDOS(bin) - 1/sqrt(bin));
    }

    overlap /= (mainWindow.nbins() - 2);

    CHECK_CLOSE(overlap, 0, 0.1);
}

TEST(RNG_MEAN)
{

    HarmOsc ho(10, 1, 1);

    double mean = 0;

    uint NC = 1E6;

    for (uint i = 0; i < NC; ++i)
    {
        mean += ho.URNG();
    }

    mean /= NC;

    CHECK_CLOSE(0.5, mean, 0.001);
}

TEST(dest)
{
    uint N = 1000;
    uint Np = 1;
    double xmax = 1.0;

    HarmOsc ho(N, Np, xmax);
    ho.changePosition(0, N/2, 1, 1);

    uint xd, d;
    for (uint dest = 0; dest < ho.freeVolume(); ++dest)
    {
        ho.findDestination(dest, xd, d, d);
        if (dest < N/2)
        {
            CHECK_EQUAL(dest, xd);
        }
        else
        {
            CHECK_EQUAL(dest + 1, xd);
        }
    }

    HarmOsc ho2(6, 4, 1);
    ho2.changePosition(3, 4, 1, 1);
    ho2.findDestination(1, xd, d, d);
    CHECK_EQUAL(5, xd);
}

TEST(RNG)
{
    uint N = 100;
    uint Np = 5;
    double xmax = 1.0;

    HarmOsc ho(N, Np, xmax);

    vec H(Np, fill::zeros);
    vec HD(N - Np, fill::zeros);
    vec HX(N, fill::zeros);

    uint NC = 1E6;

    uint pi, xd, dummy;
    for (uint i = 0; i < NC; ++i)
    {
        ho.getRandomParticleAndDestination(pi, xd, dummy, dummy);
        H(pi)++;
        uint destination = ho.URNG()*ho.freeVolume();
        HD(destination)++;
        HX(xd)++;

        ho.changePosition(pi, xd, 1, 1);
    }

    for (uint i = 0; i < Np; ++i)
    {
        if (H(i) == 0)
        {
            cout << "zero visits at H " << i << endl;
        }
    }

    for (uint i = 0; i < N; ++i)
    {
        if (N - i > Np)
        {

            if (HD(i) == 0)
            {
                cout << "zero visits at HD " << i << endl;
            }

        }

        if (HX(i) == 0)
        {
            cout << "zero visits at HX " << i << endl;
        }
    }



    CHECK_CLOSE(1.0, H.min()/mean(H), 0.1);
    CHECK_CLOSE(1.0, HD.min()/mean(HD), 0.1);
    CHECK_CLOSE(1.0, HX.min()/mean(HX), 0.1);

}

TEST(bins)
{

    using WLMC::Window;

    double min = 8;
    double max = 1024;

    uint nbins = 1024;

    HarmOsc ho(10, 1, 10);
    cout << ho.URNG() << endl;
    Window w(&ho, nbins, min, max, false);

    vec counts(nbins, fill::zeros);
    vec counts2(nbins, fill::zeros);

    double dv = 0.015625;
    vec binEnergies = linspace(min, max, nbins + 1);

    uint illegalCounter = 0;
    uint wrongBinCounter = 0;
    double deviate = 0;

    for (double value = min; value <= max; value += dv)
    {

        if (!w.isLegal(value))
        {
            illegalCounter++;
        }

        uint bfbin = 0;

        if (value == max)
        {
            bfbin = nbins - 1;
        }

        else
        {
            for (uint i = 0; i < nbins; ++i)
            {
                if (binEnergies(i) <= value && binEnergies(i + 1) > value)
                {
                    bfbin = i;
                    break;
                }
            }
        }

        uint bin = w.getBin(value);

        if (bin != bfbin)
        {

            wrongBinCounter++;
            deviate += fabs(bin - bfbin);
            cout << bin << " is not " << bfbin << " for value " << setprecision(16) << fixed << value
                 << " " << setprecision(16) << binEnergies(bfbin) << " - " << setprecision(16) << binEnergies(bfbin + 1)
                 << " | fail?: " <<  " " << setprecision(16) << binEnergies(bin) << " " << setprecision(16) << binEnergies(bin + 1) << endl;

        }
        counts(bin)++;

    }

    uint NC = 1E7;
    for (uint i = 0; i < NC; ++i)
    {
        double value = min + (max - min)*ho.URNG();

        uint bin = w.getBin(value);

        counts2(bin)++;
    }

    CHECK_CLOSE(1, counts2.min()/mean(counts2), 0.05);

    uvec indx = find(counts == 0);

    CHECK_EQUAL(0, indx.size());

    counts -= counts.min();

    //    cout << sum(counts) << endl;
    //    cout << int(sum(counts))%nbins << endl;
    //    cout << find(counts == 2) << endl;
    //    cout << find(counts == 0) << endl;
    //    cout << size(find(counts == 1)) << endl;

    if (wrongBinCounter != 0) deviate /= wrongBinCounter;

    double _mean = mean(counts);
    uint meanDeviateCounter = 0;
    for (uint bin = 0; bin < nbins; ++bin)
    {
        if (fabs(_mean - counts(bin)) > 1E-3)
        {
            meanDeviateCounter++;
        }
    }

    CHECK_EQUAL(0, illegalCounter);
    CHECK_EQUAL(0, meanDeviateCounter);
    CHECK_EQUAL(0, wrongBinCounter);
    CHECK_CLOSE(0, arma::max(arma::abs(counts - _mean)), 1E-3);
    CHECK_CLOSE(0, deviate, 1E-3);


}



int main()
{
    srand(time(NULL));

    return UnitTest::RunAllTests();
}

const WLMC::Window &harmosc()
{

    uint Np = 1;
    uint N = 1000;

    uint nbins = 40;

    double xmax = 1;

    HarmOsc ho(N, Np, xmax);

    stringstream s;
    s << Np;

    DCViz viz("stateDensity" + s.str() + ".arma");
//    viz.launch(true, 0.1, 30, 30);

    return ho.execute(nbins, 0, 1, 1E-4);

}

//void ising()
//{

//    uint nbins = 1000;

//    uint NX = 30;
//    uint NY = 30;

//    uint Np = NX*NY/2;

//    uint overlap = 20;
//    uint minWindowSize = 100;

//    ising2D system(NX, NY, Np, overlap, minWindowSize);

//    WLMC::Window *mainWindow = system.execute(nbins, 0, datum::e, 1 + 1E-6);

//    vec E = mainWindow->energies();
//    vec DOS = mainWindow->logDOS();

//    double T, dT;
//    uint N = 100000;

//    vec temperatures = linspace(0.1, 40.0, N);
//    dT = temperatures(1) - temperatures(0);

//    vec avgEs(N);

//    double partitionFunc, avgE, pdf;
//    for (uint i = 0; i < N; ++i)
//    {
//        T = temperatures(i);

//        partitionFunc = 0;
//        avgE = 0;

//        for (uint bin = 0; bin < nbins; ++bin)
//        {
//            if (mainWindow->isDeflatedBin(bin))
//            {
//                continue;
//            }

//            pdf = DOS(bin)*exp(-E(bin)/T);

//            partitionFunc += pdf;
//            avgE += pdf*E(bin);

//        }

//        avgE /= partitionFunc;

//        avgEs(i) = avgE;

//    }

//    vec C(N, fill::zeros);
//    mat CT(N, 2);

//    CT.col(1) = temperatures;


//    for (uint i = 1; i < N - 1; ++i)
//    {
//        C(i) = (avgEs(i + 1) - avgEs(i - 1))/(2*dT);
//    }

//    CT.col(0) = C;

//    CT.save("CT.arma");

//    cout << C.max()/mean(C)
//         << " at "
//         << temperatures(find(C == C.max()).eval()(0))
//         << endl;
//}
