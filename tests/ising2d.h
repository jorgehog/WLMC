#pragma once

#include <WLMC.h>
#include <BADAss/badass.h>

class ising2D : public WLMC::System
{

public:

    ising2D(const uint NX, const uint NY, const uint Np, const uint overlap, const uint minWindowSize) :
        WLMC::System(Np, NX, NY, 1, 2500, 0.85, overlap, minWindowSize, ".", [] () {return randu(1, 1).eval()(0, 0);})
    {
        m_lattice.set_size(NX, NY);
        m_lattice.zeros();

        m_particles.set_size(Np, 2);

        uint n = 0;
        for (uint x = 0; x < NX; ++x)
        {
            for (uint y = 0; y < NY; ++y)
            {
                m_particles(n, 0) = x;
                m_particles(n, 1) = y;

                m_lattice(x, y) = 1;

                n++;

                if (n == Np)
                {
                    break;
                }

            }

            if (n == Np)
            {
                break;
            }


        }

    }

    double energy(int x, int y) const
    {
        int NX = m_lattice.n_rows;
        int NY = m_lattice.n_cols;

        double E = 0;

        if (isOccupiedLoction((x + 1)%NX, y))
        {
            E++;
        }

        if (isOccupiedLoction((x - 1 + NX)%NX, y))
        {
            E++;
        }

        if (isOccupiedLoction(x, (y + 1)%NY))
        {
            E++;
        }

        if (isOccupiedLoction(x, (y - 1 + NY)%NY))
        {
            E++;
        }

        return E;

    }

    double energy(uint particleIndex) const
    {
        int x = m_particles(particleIndex, 0);
        int y = m_particles(particleIndex, 1);

        return energy(x, y);
    }

    void setMinimumConfiguration();
    void setMaximumConfiguration();
    bool isOccupiedLoction(const uint x, const uint y, const uint z = 1) const;
    double getValue(const uint particleIndex) const;
    double getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd) const;
    void changePosition(const uint particleIndex, const uint xd, const uint yd, const uint zd);
    void getPosition(const uint particleIndex, uint &x, uint &y, uint &z) const;

private:

    umat m_lattice;
    umat m_particles;
};

void ising2D::setMinimumConfiguration()
{
    uint n = 0;

    for (uint x = 0; x < m_lattice.n_rows; ++x)
    {
        for (uint y = 0; y < m_lattice.n_cols; ++y)
        {
            if (n == nParticles())
            {
                return;
            }

            if ((x + y)%2 == 0)
            {
                m_presetMinimum(n, 0) = x;
                m_presetMinimum(n, 1) = y;
                n++;
            }

        }
    }
}

void ising2D::setMaximumConfiguration()
{
    uint n = 0;

    for (uint x = 0; x < m_lattice.n_rows; ++x)
    {
        for (uint y = 0; y < m_lattice.n_cols; ++y)
        {
            if (n == nParticles())
            {
                return;
            }

            m_presetMaximum(n, 0) = x;
            m_presetMaximum(n, 1) = y;
            n++;
        }
    }
}

bool ising2D::isOccupiedLoction(const uint x, const uint y, const uint z) const
{
    (void) z;

    return m_lattice(x, y) == 1;
}

double ising2D::getValue(const uint particleIndex) const
{
    return energy(particleIndex);
}

double ising2D::getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd) const
{
    (void) zd;

    int x = m_particles(particleIndex, 0);
    int y = m_particles(particleIndex, 1);

    double diff = energy(xd, yd) - energy(x, y);

    if (((abs(signed(xd) - x) == 1 && signed(yd) == y) || (abs(signed(yd) - y) == 1 && signed(xd) == x)))
    {
        diff--;
    }

    return diff;
}

void ising2D::changePosition(const uint particleIndex, const uint xd, const uint yd, const uint zd)
{
    (void) zd;

    int x = m_particles(particleIndex, 0);
    int y = m_particles(particleIndex, 1);

    m_lattice(x, y) = 0;
    m_lattice(xd, yd) = 1;

    m_particles(particleIndex, 0) = xd;
    m_particles(particleIndex, 1) = yd;

}

void ising2D::getPosition(const uint particleIndex, uint &x, uint &y, uint &z) const
{
    x = m_particles(particleIndex, 0);
    y = m_particles(particleIndex, 1);
    z = 0;

}
