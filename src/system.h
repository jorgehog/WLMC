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
           function<double()> URNG);

    virtual ~System() {}

    virtual void setMaximumConfiguration();

    virtual void setMinimumConfiguration();

    virtual bool isOccupiedLoction(const uint x, const uint y, const uint z) const = 0;

    virtual double getValue(const uint particleIndex) const = 0;

    virtual double getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd);

    virtual double getTotalValue() const;

    virtual void changePosition(const uint particleIndex, const uint xd, const uint yd, const uint zd) = 0;

    virtual void getPosition(const uint particleIndex, uint &x, uint &y, uint &z) const = 0;

    virtual void savePositionData(const uint framenumber) const;

    void consistencyCheckOptimizedValues();

    Window *execute(const uint nbins, const double adaptive, const double fStart, const double fEnd, function<double(double)> reduceFunction);

    Window *execute(const uint nbins, const double adaptive, const double fStart, const double fEnd)
    {
        return execute(nbins, adaptive, fStart, fEnd, [] (double pre) {return sqrt(pre);});
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

    void clipWindow(Window &window) const;

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

    const double &f() const
    {
        return m_f;
    }

    const string path() const
    {
        return m_path;
    }

    enum class extrema
    {
        minimum,
        maximum
    };

private:

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

    double m_f;

    const function<double()> m_URNG;

    const string m_path;

    ucube m_presetWindowConfigurations;
    vec m_presetWindowValues;

    double getGlobalExtremum(const extrema type);

    void getRandomParticleAndDestination(uint &particleIndex, uint &xd, uint &yd, uint &zd) const;

    void randomizeParticlePositions();

protected:

    umat m_presetMinimum;
    umat m_presetMaximum;


};

}
