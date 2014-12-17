#pragma once

#include <armadillo>

#include <functional>

using namespace arma;
using namespace std;

namespace WLMC
{

class Window;

class System
{
public:
    System(const uint nParticles,
           const uint NX,
           const uint NY,
           const uint NZ,
           const uint movesPerSampling,
           const double flatnessCriterion,
           const uint overlap,
           const uint minWindowSize,
           const double flatnessGradientTreshold,
           const double deflationLimit,
           const string path,
           function<double()> URNG,
           bool requiresDestination = true,
           bool requireParticlePositions = true);

    virtual ~System();

    virtual void onSuggestedTrialMove(const uint particleIndex, const uint xd, const uint yd, const uint zd)
    {
        (void) particleIndex;

        (void) xd;
        (void) yd;
        (void) zd;
    }

    virtual void onPresetSave(const uint n)
    {
        (void) n;
    }

    virtual void onPresetLoad(const uint n)
    {
        (void) n;
    }

    uint numberOfPresets(const uint nbins) const;

    virtual double setMaximumConfiguration();

    virtual double setMinimumConfiguration();

    virtual bool isOccupiedLoction(const uint x, const uint y, const uint z) const = 0;

    virtual double getValue(const uint particleIndex) const = 0;

    virtual double getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd);

    virtual double getTotalValue() const;

    virtual void changePosition(const uint particleIndex, const uint xd, const uint yd, const uint zd) = 0;

    virtual void getPosition(const uint particleIndex, uint &x, uint &y, uint &z) const = 0;

    virtual void savePositionData(const uint framenumber) const;

    virtual void shiftDOS(vec &logDOS)
    {
        double m = logDOS.max();
        logDOS -= m;
    }

    void consistencyCheckOptimizedValues();

    Window &execute(const uint nbins, const double adaptive, const double fStart, const double fEnd, function<double(double)> reduceFunction);

    Window &execute(const uint nbins, const double adaptive, const double logfStart, const double logfEnd, const double logfReduceFactor = 0.5)
    {
        return execute(nbins, adaptive, logfStart, logfEnd, [&logfReduceFactor] (double pre) {return logfReduceFactor*pre;});
    }

    void sampleWindow(Window *window);

    bool doWLMCMove(Window *window);

    void doRandomMove();

    void findDestination(const uint destination, uint &xd, uint &yd, uint &zd) const;

    void locateGlobalExtremaValues(double &min, double &max);

    void setupPresetWindowConfigurations(Window &mainWindow);

    void loadConfigurationForWindow(const Window *window);

    void loadConfiguration(const umat &config);

    uint getPresetBinFromValue(const double value) const;

    void getRandomParticleAndDestination(uint &particleIndex, uint &xd, uint &yd, uint &zd) const;

    void clipWindow(Window &window) const;

    double URNG() const
    {
        return m_URNG();
    }

    const uint &freeVolume() const
    {
        return m_freeVolume;
    }

    const uint &nParticles() const
    {
        return m_nParticles;
    }

    const uint &movesPerWindowCheck() const
    {
        return m_movesPerSampling;
    }

    const double &flatnessCriterion() const
    {
        return m_flatnessCriterion;
    }

    const double &flatnessGradientTreshold() const
    {
        return m_flatnessGradientTreshold;
    }

    const double &deflationLimit() const
    {
        return m_deflationLimit;
    }

    const uint &overlap() const
    {
        return m_overlap;
    }

    uint minWindowSize() const
    {
        return m_minWindowSize;
    }

    const double &logf() const
    {
        return m_logf;
    }

    const string path() const
    {
        return m_path;
    }

    const uint &NX() const
    {
        return m_NX;
    }

    const uint &NY() const
    {
        return m_NY;
    }

    const uint &NZ() const
    {
        return m_NZ;
    }

    enum class extrema
    {
        minimum,
        maximum
    };

private:

    Window *m_mainWindow;

    const uint m_nParticles;

    const uint m_NX;
    const uint m_NY;
    const uint m_NZ;

    const uint m_volume;
    const uint m_freeVolume;

    const uint m_movesPerSampling;

    const double m_flatnessCriterion;
    const double m_deflationLimit;
    const uint m_overlap;
    const uint m_minWindowSize;

    const double m_flatnessGradientTreshold;

    double m_logf;

    const function<double()> m_URNG;

    const string m_path;

    const bool m_requiresDestination;
    const bool m_requiresPositions;

    ucube m_presetWindowConfigurations;
    vec m_presetWindowValues;

    double getGlobalExtremum(const extrema type);

    void randomizeParticlePositions();

    umat m_presetMinimum;
    umat m_presetMaximum;


};

}
