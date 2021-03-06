#pragma once

#include <armadillo>
#include <vector>

using std::vector;

using namespace arma;

namespace WLMC
{

class System;

struct WindowParams;

class Window
{
public:

    enum class OverlapTypes
    {
        Lower,
        Upper,
        None
    };

    Window(System *system,
           const vec &parentLogDOS,
           const vec &parentEnergies, const vector<bool> &parentBinDeflation,
           const uint lowerLimitOnParent,
           const uint upperLimitOnParent,
           Window::OverlapTypes overlapType,
           bool allowSubwindowing);

    Window(System *system,
           const uint nBins,
           const double minValue,
           const double maxValue,
           bool allowSubwindowing);

    Window(const Window &parentWindow, const WindowParams &windowParams);

    virtual ~Window();

    void adapt(const uint lowerLimit, const uint upperLimit);

    vec getHistogram(const uint nmoves);

    void loadInitialConfig();

    void calculateWindow();

    double estimateFlatness(const uint lowerLimit, const uint upperLimit) const;

    double estimateFlatness() const
    {
        return estimateFlatness(0, m_nbins);
    }

    void findSubWindows();

    double getMeanFlatness(const uint lowerLimit, const uint upperLimit) const;

    double getMeanFlatness() const
    {
        return getMeanFlatness(0, m_nbins);
    }

    void findComplementaryRoughAreas(vector<WindowParams> &roughWindowParams) const;

    bool findFlatArea();

    bool scanForFlattestArea();

    double sparsity(const uint lowerLimit, const uint upperLimit) const;

    void expandFlattestArea();

    bool flatProfileIsContinousOnParent() const;


    void getSubWindowLimits(WindowParams &windowParams) const;

    void registerVisit(const uint bin);

    uint getBin(double value) const;

    void reset();

    void resetDOS();

    void deflateBin(const uint bin)
    {
        cout << "deflated " << bin << " " << logDOS(bin) << endl;
        m_visitCounts(bin) = m_unsetCount;
        m_deflatedBins.at(bin) = true;
    }

    bool allowsSubwindowing() const;

    bool flatspanGradientConverged() const;

    bool isLegal(const double value) const
    {
        return value >= m_minValue && value <= m_maxValue;
    }

    const double &minValue() const
    {
        return m_minValue;
    }

    const double &maxValue() const
    {
        return m_maxValue;
    }

    const double &valueSpan() const
    {
        return m_valueSpan;
    }

    const uint &lowerLimitOnParent() const
    {
        return m_lowerLimitOnParent;
    }

    void rofl(const uint p)
    {
        m_lowerLimitOnParent = m_upperLimitOnParent = p;
    }

    const uint &upperLimitOnParent() const
    {
        return m_upperLimitOnParent;
    }

    System *system() const
    {
        return m_system;
    }

    const vec &logDOS() const
    {
        return m_logDOS;
    }

    const double &logDOS(const uint i) const
    {
        return m_logDOS(i);
    }

    void logDOS(const vec newLogDOS)
    {
        m_logDOS = newLogDOS;
    }

    const vec &energies() const
    {
        return m_energies;
    }

    const uvec &visitCounts() const
    {
        return m_visitCounts;
    }

    const uint &visitCounts(const uint i) const
    {
        return m_visitCounts(i);
    }

    const vector<bool> &deflatedBins() const
    {
        return m_deflatedBins;
    }

    bool isDeflatedBin(const uint bin) const
    {
        return m_deflatedBins.at(bin);
    }

    const uint &nbins() const
    {
        return m_nbins;
    }

    bool isUnsetCount(const uint i) const
    {
        return m_visitCounts(i) == m_unsetCount;
    }

    bool isFlat(const uint lowerLimit, const uint upperLimit) const;

    bool isFlat() const
    {
        return isFlat(0, m_nbins);
    }

    bool isFlatOnParent() const;

    const Window::OverlapTypes &overlapType() const
    {
        return m_overlapType;
    }

    bool overlapsAtTop() const
    {
        return m_overlapType == OverlapTypes::Upper;
    }

    bool overlapsAtBottom() const
    {
        return m_overlapType == OverlapTypes::Lower;
    }

    void dump_output() const;

    bool atBoundaryValue() const;

    bool atBoundaryValue(double &boundaryValue) const;

    static constexpr uint m_unsetCount = std::numeric_limits<uint>::max();

private:

    System *m_system;

    vector<Window*> m_subWindows;
    bool m_allowsSubwindowing;

    const Window::OverlapTypes m_overlapType;

    uint m_lowerLimitOnParent;
    uint m_upperLimitOnParent;
    uint m_nbins;

    double m_minValue;
    double m_maxValue;
    double m_valueSpan;

    vec m_logDOS;
    vec m_energies;
    uvec m_visitCounts;

    uint m_flatAreaLower;
    uint m_flatAreaUpper;

    uint m_gradientSampleCounter;

    double m_spanSum;
    vec4 m_spanSums;
    double m_spanGradient;
    double m_spanLaplace;    

    double m_centerSum;
    vec4 m_centerSums;
    double m_centerGradient;
    vector<bool> m_deflatedBins;


    uint m_outputLevel;

    void deflateDOS();

    void normaliseDOS()
    {
        double m = m_logDOS.max();
        m_logDOS -= m;

        deflateDOS();
    }

    void mergeWith(Window *other);

    uint getOverlapPoint(const Window *other);

};

}
