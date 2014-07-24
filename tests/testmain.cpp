#include <WLMC.h>
#include <armadillo>
#include <vector>
#include <algorithm>

using namespace arma;
using namespace std;

class ising2D : public WLMC::System
{

public:

    ising2D(const uint NX, const uint NY, const uint Np, const uint overlap, const uint minWindowSize) :
        WLMC::System(Np, NX, NY, 1, 1000, 0.85, overlap, minWindowSize, ".", [] () {return randu(1, 1).eval()(0, 0);})
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

    bool isOccupiedLoction(const uint x, const uint y, const uint z = 1) const;
    double getValue(const uint particleIndex) const;
    double getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd) const;
    double getTotalValue() const;
    void changePosition(const uint particleIndex, const uint xd, const uint yd, const uint zd);
    void getPosition(const uint particleIndex, uint &x, uint &y, uint &z) const;

private:

    umat m_lattice;
    umat m_particles;
};

bool ising2D::isOccupiedLoction(const uint x, const uint y, const uint z) const
{
    (void) z;

    return m_lattice(x, y) == 1;
}

double ising2D::getValue(const uint particleIndex) const
{
    int x = m_particles(particleIndex, 0);
    int y = m_particles(particleIndex, 1);

    return energy(x, y);
}

double ising2D::getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd) const
{
    (void) zd;

    int x = m_particles(particleIndex, 0);
    int y = m_particles(particleIndex, 1);

    return energy(xd, yd) - energy(x, y);
}

double ising2D::getTotalValue() const
{
    double E = 0;

    for (uint p = 0; p < m_particles.n_rows; ++p)
    {
        int x = m_particles(p, 0);
        int y = m_particles(p, 1);

        E += energy(x, y);
    }

    return E;
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
    z = 1;

}


int main()
{
    uint nbins = 100;

    uint NX = 10;
    uint NY = 10;

    uint Np = 20;

    uint overlap = nbins/20;
    uint minWindowSize = nbins/10;

    ising2D system(NX, NY, Np, overlap, minWindowSize);

    system.execute(nbins, 0, datum::e, 1 + 1E-6);

    return 0;
}
