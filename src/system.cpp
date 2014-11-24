#include "system.h"
#include "window.h"
#include "windowparams.h"

#include <iomanip>

#include <BADAss/badass.h>
#include <lammpswriter/lammpswriter.h>

using namespace WLMC;

System::System(const uint nParticles,
               const uint NX,
               const uint NY,
               const uint NZ,
               const uint movesPerSampling,
               const double flatnessCriterion,
               const uint overlap,
               const uint minWindowSizeRough,
               const double flatnessGradientTreshold,
               const double deflationLimit,
               const string path,
               function<double()> URNG,
               bool requiresDestination) :
    m_nParticles(nParticles),
    m_NX(NX),
    m_NY(NY),
    m_NZ(NZ),
    m_volume(m_NX*m_NY*m_NZ),
    m_freeVolume(m_volume - nParticles),
    m_movesPerSampling(movesPerSampling),
    m_flatnessCriterion(flatnessCriterion),
    m_deflationLimit(deflationLimit),
    m_overlap(overlap),
    m_minWindowSize(minWindowSizeRough),
    m_flatnessGradientTreshold(flatnessGradientTreshold),
    m_URNG(URNG),
    m_path(path),
    m_requiresDestination(requiresDestination)
{
    m_presetMinimum.set_size(m_nParticles, 3);
    m_presetMaximum.set_size(m_nParticles, 3);

    BADAss(nParticles, >=, 1, "Number of particles.");
    BADAss(nParticles, <=, m_volume, "Incorrect number of particles.");

}

System::~System()
{
    delete m_mainWindow;
}

uint System::numberOfPresets(const uint nbins) const
{
    return ceil(nbins/double(m_minWindowSize + m_overlap));
}

double System::setMaximumConfiguration()
{
    uint sweep, nSweeps, x, y, z;;
    double localMax;
    bool isIn;

    vector<double> allExtrema;

    double max = std::numeric_limits<double>::min();

    sweep = 0;
    nSweeps = 10;

    while (sweep < nSweeps)
    {
        cout << "sweep " << sweep + 1 << " / " << nSweeps << endl;

        randomizeParticlePositions();

        localMax = getGlobalExtremum(System::extrema::maximum);

        isIn = false;
        for (double extrema : allExtrema)
        {
            if (fabs(localMax - extrema) < 0.001)
            {
                isIn = true;
                break;
            }
        }
        if (!isIn)
        {
            allExtrema.push_back(localMax);
        }

        if (localMax > max)
        {
            max = localMax;

            for (uint particleIndex = 0; particleIndex < m_nParticles; ++particleIndex)
            {
                getPosition(particleIndex, x, y, z);

                m_presetMaximum(particleIndex, 0) = x;
                m_presetMaximum(particleIndex, 1) = y;
                m_presetMaximum(particleIndex, 2) = z;
            }
        }

        sweep++;
    }

    loadConfiguration(m_presetMaximum);

    return max;

}

double System::setMinimumConfiguration()
{
    uint nSweeps, sweep, x, y, z;
    double localMin;
    bool isIn;
    vector<double> allExtrema;

    double min = std::numeric_limits<double>::max();

    nSweeps = 10;
    sweep = 0;

    while (sweep < nSweeps)
    {
        cout << "sweep " << sweep + 1 << " / " << nSweeps << endl;

        randomizeParticlePositions();

        localMin = getGlobalExtremum(System::extrema::minimum);

        isIn = false;
        for (double extrema : allExtrema)
        {
            if (fabs(localMin - extrema) < 0.001)
            {
                isIn = true;
                break;
            }
        }
        if (!isIn)
        {
            allExtrema.push_back(localMin);
        }

        if (localMin < min)
        {
            min = localMin;

            for (uint particleIndex = 0; particleIndex < m_nParticles; ++particleIndex)
            {
                getPosition(particleIndex, x, y, z);

                m_presetMinimum(particleIndex, 0) = x;
                m_presetMinimum(particleIndex, 1) = y;
                m_presetMinimum(particleIndex, 2) = z;
            }
        }

        sweep++;
    }

    loadConfiguration(m_presetMinimum);

    return min;

}

double System::getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd)
{
    uint xPre, yPre, zPre;
    getPosition(particleIndex, xPre, yPre, zPre);

    changePosition(particleIndex, xd, yd, zd);

    double valDest = getTotalValue();

    changePosition(particleIndex, xPre, yPre, zPre);

    return valDest - getTotalValue();

}

double System::getTotalValue() const
{
    double totalValue = 0;

    for (uint particleIndex = 0; particleIndex < m_nParticles; ++particleIndex)
    {
        totalValue += getValue(particleIndex);
    }

    return totalValue;
}

void System::savePositionData(const uint framenumber) const
{
    uint x, y, z;

    lammpswriter writer(5, "WLMC", m_path);

    writer.setSystemSize(m_NX, m_NY, m_NZ);

    writer.initializeNewFile(framenumber, m_nParticles);

    for (uint particleIndex = 0; particleIndex < m_nParticles; ++particleIndex)
    {
        getPosition(particleIndex, x, y, z);

        writer << 0 << x << y << z << getValue(particleIndex);
    }

    writer.finalize();
}

void System::consistencyCheckOptimizedValues()
{
#ifndef NDEBUG

    uint particleIndex, x, y, z, xd, yd, zd;

    double value = getTotalValue();

    double valueCumSum = 0;
    for (particleIndex = 0; particleIndex < nParticles(); ++particleIndex)
    {
        valueCumSum += getValue(particleIndex);
    }

    BADAssClose(value, valueCumSum, 1E-10, "TotalValue and local values mismatch when accumulated.");

    getRandomParticleAndDestination(particleIndex, xd, yd, zd);


    double valueDifference = getValueDifference(particleIndex, xd, yd, zd);

    getPosition(particleIndex, x, y, z);

    BADAssBool(isOccupiedLoction(x, y, z), "Error in occupancy check.");
    BADAssBool(!isOccupiedLoction(xd, yd, zd), "Error in occupancy check.");

    changePosition(particleIndex, xd, yd, zd);

    BADAssBool(!isOccupiedLoction(x, y, z), "Error in occupancy check.");
    BADAssBool(isOccupiedLoction(xd, yd, zd), "Error in occupancy check.");

    double valueDifferenceBruteForce = getTotalValue() - value;

    changePosition(particleIndex, x, y, z);

    BADAssClose(valueDifferenceBruteForce, valueDifference, 1E-10, "valueDifference implementation mismatches brute force calculation.");

#endif
}

const Window &System::execute(const uint nbins, const double adaptive, const double logfStart, const double logfEnd, function<double(double)> reduceFunction)
{
//    consistencyCheckOptimizedValues();

    m_presetWindowConfigurations.set_size(m_nParticles, 3, numberOfPresets(nbins));

    double max, min;
    locateGlobalExtremaValues(min, max);

    m_mainWindow = new Window(this, nbins, min, max, adaptive);

    setupPresetWindowConfigurations(*m_mainWindow);

    m_logf = logfStart;
    while (m_logf >= logfEnd)
    {
        m_mainWindow->calculateWindow();

        m_logf = reduceFunction(m_logf);

        m_mainWindow->reset();
    }

    m_mainWindow->dump_output();
    return *m_mainWindow;

}

void System::sampleWindow(Window *window)
{
    uint nMoves = 0;
    
    while (nMoves < m_movesPerSampling)
    {
        if (doWLMCMove(window))
        {
            nMoves++;
            BADAssBool(window->isLegal(getTotalValue()));
        }
    }
}

bool System::doWLMCMove(Window *window)
{
    uint particleIndex, xd, yd, zd;

    getRandomParticleAndDestination(particleIndex, xd, yd, zd);

    double oldValue = getTotalValue();
    uint oldBin = window->getBin(oldValue);

    onSuggestedTrialMove(particleIndex, xd, yd, zd);

    double newValue = oldValue + getValueDifference(particleIndex, xd, yd, zd);

    if (!window->isLegal(newValue))
    {
        return false;
    }

    uint newBin = window->getBin(newValue);

    if (window->isDeflatedBin(oldBin) || window->isDeflatedBin(newBin))
    {
        changePosition(particleIndex, xd, yd, zd);
        return false;
    }


    double oldLogDOS = window->logDOS(oldBin);
    double newLogDOS = window->logDOS(newBin);

    bool accepted = true;

    if (oldLogDOS < newLogDOS)
    {
        accepted = (m_URNG() < exp(oldLogDOS - newLogDOS));
    }

    if (accepted)
    {
        changePosition(particleIndex, xd, yd, zd);

        window->registerVisit(newBin);
    }

    else
    {
        window->registerVisit(oldBin);
    }

    return true;

}

void System::doRandomMove()
{
    uint particleIndex, xd, yd, zd;

    xd = yd = zd = 0;

    getRandomParticleAndDestination(particleIndex, xd, yd, zd);

    onSuggestedTrialMove(particleIndex, xd, yd, zd);

    changePosition(particleIndex, xd, yd, zd);
}

void System::findDestination(const uint destination, uint &xd, uint &yd, uint &zd) const
{
    uint search = 0;

    xd = yd = zd = 0;

    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                if (!isOccupiedLoction(x, y, z))
                {
                    if (search == destination)
                    {
                        xd = x;
                        yd = y;
                        zd = z;

                        return;
                    }

                    search++;
                }
            }
        }
    }

    BADAssBreak("Unable to find destination", [&] ()
    {
        cout << "dest = " << destination << " freevolume = " << m_freeVolume << endl;
    });

}

void System::locateGlobalExtremaValues(double &min, double &max)
{
    uint x, y, z, n;

    n = m_presetWindowConfigurations.n_slices;

    min = setMinimumConfiguration();

    for (uint i = 0; i < nParticles(); ++i)
    {
        getPosition(i, x, y, z);
        m_presetWindowConfigurations(i, 0, 0) = x;
        m_presetWindowConfigurations(i, 1, 0) = y;
        m_presetWindowConfigurations(i, 2, 0) = z;
    }

    onPresetSave(0);

    max = setMaximumConfiguration();

    for (uint i = 0; i < nParticles(); ++i)
    {
        getPosition(i, x, y, z);
        m_presetWindowConfigurations(i, 0, n - 1) = x;
        m_presetWindowConfigurations(i, 1, n - 1) = y;
        m_presetWindowConfigurations(i, 2, n - 1) = z;
    }

    onPresetSave(n - 1);

    cout << "located extrema " << min << " " << max << endl;
}

void System::setupPresetWindowConfigurations(Window &mainWindow)
{
    if (!mainWindow.allowsSubwindowing())
    {
        return;
    }

    uint x, y, z, bin, nSet;
    double value, max, min;

    max = mainWindow.maxValue();
    min = mainWindow.minValue();

    uint n = numberOfPresets(mainWindow.nbins());

    m_presetWindowValues = linspace(min, max, n + 1);

    uvec binSet = zeros<uvec>(n);
    binSet(0) = 1;
    binSet(n - 1) = 1;
    nSet = 2;

    uint nMax = 10000;
    uint nMax2 = 1000;
    uint nTest = 0;

    double spaciousness = 10;

    while (nSet < n)
    {
        doRandomMove();

        nTest++;
        if (nTest > nMax)
        {
            for (uint bin = 0; bin < n; ++bin)
            {
                if (binSet(bin))
                {
                    continue;
                }

                uint lower = mainWindow.getBin(m_presetWindowValues(bin));
                uint upper = mainWindow.getBin(m_presetWindowValues(bin + 1));

                cout << "unable to find configuration in window "
                     << m_presetWindowValues(bin) << " "
                     << m_presetWindowValues(bin + 1)
                     << ". Deflating." << endl;

                for (uint mwBin = lower; mwBin < upper; ++mwBin)
                {
                    mainWindow.deflateBin(mwBin);
                }
            }

            loadConfiguration(m_presetWindowConfigurations.slice(0));
            onPresetLoad(0);
            savePositionData(4);

            return;
        }

        value = getTotalValue();

        //Just in case the extrema calculations went bad
        if (value < min || value > max)
        {
            if ((nTest + 1)%nMax2 == 0)
            {
                loadConfiguration(m_presetMinimum);
                onPresetLoad(0);
            }

            continue;
        }

        bin = getPresetBinFromValue(value);

        double mmin = m_presetWindowValues(bin);
        double mmax = m_presetWindowValues(bin + 1);

        double delta = (mmax - mmin)/spaciousness;

        bool valueIsTooCloseToBoundary = value < mmin + delta || value > mmax - delta;

        if (binSet(bin) != 0 || valueIsTooCloseToBoundary)
        {
            continue;
        }

        binSet(bin) = 1;
        nTest = 0;

        for (uint particleIndex = 0; particleIndex < m_nParticles; ++particleIndex)
        {
            getPosition(particleIndex, x, y, z);

            m_presetWindowConfigurations(particleIndex, 0, bin) = x;
            m_presetWindowConfigurations(particleIndex, 1, bin) = y;
            m_presetWindowConfigurations(particleIndex, 2, bin) = z;
        }

        onPresetSave(bin);

        cout << "set energy for value " << value << " in bin " << bin << " with energyspan " << m_presetWindowValues(bin) << " - " << m_presetWindowValues(bin + 1) << " " << nSet+1 << "/" << n << endl;

        nSet++;

        loadConfiguration(m_presetMinimum);
        onPresetLoad(0);

    }

}

void System::loadConfigurationForWindow(const Window *window)
{

    double value = (window->minValue() + window->maxValue())/2;

    uint bin = getPresetBinFromValue(value);

    cout << "value " << value << " belongs in bin " << bin << " between " << m_presetWindowValues(bin) << " - " << m_presetWindowValues(bin + 1) << endl;

    loadConfiguration(m_presetWindowConfigurations.slice(bin));
    onPresetLoad(bin);

}

void System::loadConfiguration(const umat &config)
{
    uint xPreset, yPreset, zPreset, x, y, z, xAvailable, yAvailable, zAvailable, particleIndexConfig;

    xAvailable = yAvailable = zAvailable = 0;

    bool isAlreadyOccupied;

    for (uint particleIndex = 0; particleIndex < m_nParticles; ++particleIndex)
    {
        getPosition(particleIndex, x, y, z);

        isAlreadyOccupied = false;
        for (particleIndexConfig = 0; particleIndexConfig < m_nParticles; ++particleIndexConfig)
        {
            xPreset = config(particleIndexConfig, 0);
            yPreset = config(particleIndexConfig, 1);
            zPreset = config(particleIndexConfig, 2);

            if (xPreset == x && yPreset == y && zPreset == z)
            {
                isAlreadyOccupied = true;
                break;
            }

            else if (!isOccupiedLoction(xPreset, yPreset, zPreset))
            {
                xAvailable = xPreset;
                yAvailable = yPreset;
                zAvailable = zPreset;
            }
        }


        if (!isAlreadyOccupied)
        {
            BADAssBool(!isOccupiedLoction(xAvailable, yAvailable, zAvailable), "Error: Site should be available.", [&] ()
            {
                cout << xAvailable << " " << yAvailable << " " << zAvailable << endl;

                cout << config << endl;
            });
            changePosition(particleIndex, xAvailable, yAvailable, zAvailable);
        }
    }

}

uint System::getPresetBinFromValue(const double value) const
{
    uint bin = 0;
    while (true)
    {
        if (value >= m_presetWindowValues(bin) && value <= m_presetWindowValues(bin + 1))
        {
            return bin;
        }

        bin++;
    }

}

void System::clipWindow(Window &window) const
{

    uint histSamples, upperLimit, lowerLimit;

    if (m_movesPerSampling > 10000)
    {
        histSamples = 10000;
    }
    else
    {
        histSamples = m_movesPerSampling;
    }

    vec hist = window.getHistogram(histSamples);

    double m = arma::max(hist);

    uint clipSize = 10;
    double thresh = 1E-2;

    uint upperLimitClip = window.nbins();
    uint lowerLimitClip = window.nbins() - clipSize;

    while (mean(hist(span(lowerLimitClip, upperLimitClip - 1))) < thresh*m)
    {
        upperLimitClip -= clipSize;
        lowerLimitClip -= clipSize;
    }

    upperLimit = (lowerLimitClip + upperLimitClip)/2;

    lowerLimitClip = 0;
    upperLimitClip = clipSize;

    while (mean(hist(span(lowerLimitClip, upperLimitClip - 1))) < thresh*m)
    {
        upperLimitClip += clipSize;
        lowerLimitClip += clipSize;
    }

    lowerLimit = (lowerLimitClip + upperLimitClip)/2;

    cout << "should calc from " << lowerLimit << " to " << upperLimit << " ? " << endl;

    hist.save("/tmp/hist.arma");

    window.adapt(lowerLimit, upperLimit);

}

double System::getGlobalExtremum(const System::extrema type)
{
    using std::min_element;
    typedef std::function<bool(const double &, const double &)> compFuncType;


    double localExtrema, valueDifference;
    uint xd, yd, zd, particleIndex;

    xd = yd = zd = 0;

    //vector representing particles = 0, 1, ..., nParticles - 1
    vector<uint> particleIndices(m_nParticles);
    for (uint particleIndex = 0; particleIndex < m_nParticles; ++particleIndex)
    {
        particleIndices.at(particleIndex) = particleIndex;
    }

    compFuncType lessThan = [] (const double &v1, const double & v2) {return v1 < v2;};
    compFuncType greaterThan = [] (const double &v1, const double & v2) {return v1 > v2;};

    compFuncType arrangeSortCompare;
    compFuncType extremumCheckCompare;

    if (type == System::extrema::maximum)
    {
        extremumCheckCompare = greaterThan;
        arrangeSortCompare = lessThan;
    }
    else
    {
        extremumCheckCompare = lessThan;
        arrangeSortCompare = greaterThan;
    }

    localExtrema = 1.0;
    double ePrev;
    while (localExtrema != 0)
    {

        ePrev = getTotalValue();
        localExtrema = 0;

        //Arrange particles by values
        particleIndex = *(std::min_element(particleIndices.begin(),
                                           particleIndices.end(),
                                           [this, arrangeSortCompare] (const uint &p1, const uint &p2)
        {
            return arrangeSortCompare(getValue(p1), getValue(p2));
        }));

        //Search for the displacement which extremizes value gain
        for (uint x = 0; x < m_NX; ++x)
        {
            for (uint y = 0; y < m_NY; ++y)
            {
                for (uint z = 0; z < m_NZ; ++z)
                {
                    if (isOccupiedLoction(x, y, z))
                    {
                        continue;
                    }

                    valueDifference = getValueDifference(particleIndex, x, y, z);

                    if (extremumCheckCompare(valueDifference, localExtrema))
                    {
                        localExtrema = valueDifference;

                        xd = x;
                        yd = y;
                        zd = z;
                    }
                }
            }
        }

        //if no such displacement exist, we go on to the next particle.
        //if no such displacement exist for all particles, then we end.
        //if it does, we perform it, then reset the particle loop.
        if (localExtrema == 0)
        {
            break;
        }
        else
        {
            changePosition(particleIndex, xd, yd, zd);

            if (fabs(getTotalValue() - ePrev) < 1E-10)
            {
                break;
            }
        }

    }

    return getTotalValue();

}

void System::getRandomParticleAndDestination(uint &particleIndex, uint &xd, uint &yd, uint &zd) const
{
    if (!m_requiresDestination)
    {
        particleIndex = m_URNG()*m_volume;
        xd = yd = zd = 0;
        return;
    }

    particleIndex = m_URNG()*m_nParticles;
    uint destination = m_URNG()*m_freeVolume;

    findDestination(destination, xd, yd, zd);
}

void System::randomizeParticlePositions()
{
    uint xd, yd, zd, destination;

    xd = yd = zd = 0;

    uint nSweeps = 5;
    uint sweep = 0;

    while (sweep < nSweeps)
    {
        for (uint particleIndex = 0; particleIndex < m_nParticles; ++particleIndex)
        {

            destination = m_URNG()*m_freeVolume;

            findDestination(destination, xd, yd, zd);

            changePosition(particleIndex, xd, yd, zd);

        }

        sweep++;

    }

}

