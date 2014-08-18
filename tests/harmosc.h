#pragma once

#include <WLMC.h>
#include <BADAss/badass.h>

class HarmOsc : public WLMC::System
{
public:

    HarmOsc(const uint N, const uint Np, const double xmax) :
        WLMC::System(Np, N, 1, 1, 10000, 0.95, 5, 10, 1.0, 1E-10, ".", [] () {return randu(1, 1).eval()(0, 0);}),
        m_xmax(xmax)
    {
        m_lattice.set_size(N);
        m_lattice.zeros();

        m_nParticles.set_size(Np);

        for (uint i = 0; i < Np; ++i)
        {
            m_nParticles(i) = i;
        }

        for (const uint &x : m_nParticles)
        {
            m_lattice(x) = 1;
        }
    }

    // System interface
public:

    double xmap(const uint lx) const
    {
        return lx*m_xmax/m_lattice.size();
    }

    double energy(const uint particleIndex) const
    {
        uint lx = m_nParticles(particleIndex);
        double x = xmap(lx);

        return x*x;
    }

    bool isOccupiedLoction(const uint x, const uint y, const uint z) const
    {
        (void) y;
        (void) z;

        return m_lattice(x) == 1;

    }
    double getValue(const uint particleIndex) const
    {
        return energy(particleIndex);
    }
    void changePosition(const uint particleIndex, const uint xd, const uint yd, const uint zd)
    {
        (void) yd;
        (void) zd;

        uint x = m_nParticles(particleIndex);

        m_lattice(x) = 0;
        m_lattice(xd) = 1;

        m_nParticles(particleIndex) = xd;

    }
    void getPosition(const uint particleIndex, uint &x, uint &y, uint &z) const
    {
        x = m_nParticles(particleIndex);
        y = z = 0;
    }

    uvec m_nParticles;
    uvec m_lattice;
    double m_xmax;

    // System interface
public:
    void setMaximumConfiguration()
    {
        m_presetMaximum(0, 0) = m_lattice.size() - 1;
    }
    void setMinimumConfiguration()
    {
        m_presetMinimum(0, 0) = 0;
    }
};
