//#pragma once

//#include <WLMC.h>
//#include <BADAss/badass.h>

//class ising2D : public WLMC::System
//{

//public:

//    ising2D(const uint L, const uint Np, const uint overlap, const uint minWindowSize) :
//        WLMC::System(Np, L, L, 1, 1000, 0.85, overlap, minWindowSize, 1.0, 1E-100, ".", [] () {return randu(1, 1).eval()(0, 0);}, false)
//    {
//        m_lattice.set_size(L, L);
//        m_lattice.fill(-1);




//    }

//    double energy(int x, int y) const
//    {
//        int NX = m_lattice.n_rows;
//        int NY = m_lattice.n_cols;

//        double E = 0;

//        E += m_lattice((x + 1)%NX, y)*m_lattice(x, y);      //left
//        E += m_lattice((x - 1 + NX)%NX, y)*m_lattice(x, y); //right
//        E += m_lattice(x, (y + 1)%NY)*m_lattice(x, y);      //top
//        E += m_lattice(x, (y - 1 + NY)%NY)*m_lattice(x, y); //bottom

//        return E;

//    }

//    double energy(uint particleIndex) const
//    {
//        int x = m_particles(particleIndex, 0);
//        int y = m_particles(particleIndex, 1);

//        return energy(x, y);
//    }

//    bool isOccupiedLoction(const uint x, const uint y, const uint z = 1) const;
//    double getValue(const uint particleIndex) const;
////    double getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd) const;
//    void changePosition(const uint particleIndex, const uint xd, const uint yd, const uint zd);
//    void getPosition(const uint particleIndex, uint &x, uint &y, uint &z) const;

//private:

//    imat m_lattice;
//};

//bool ising2D::isOccupiedLoction(const uint x, const uint y, const uint z) const
//{
//    (void) x;
//    (void) y;
//    (void) z;

//    return false;
//}

//double ising2D::getValue(const uint particleIndex) const
//{
//    return energy(particleIndex);
//}

////double ising2D::getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd) const
////{
////    (void) zd;

////    int x = m_particles(particleIndex, 0);
////    int y = m_particles(particleIndex, 1);

////    double diff = energy(xd, yd) - energy(x, y);

////    if (((abs(signed(xd) - x) == 1 && signed(yd) == y) || (abs(signed(yd) - y) == 1 && signed(xd) == x)))
////    {
////        diff--;
////    }

////    return diff;
////}

//void ising2D::changePosition(const uint particleIndex, const uint xd, const uint yd, const uint zd)
//{
//    (void) xd;
//    (void) yd;
//    (void) zd;

//    int x = particleIndex;
//    int y = m_particles(particleIndex, 1);

//    m_lattice(x, y) *= -1;
//}

//void ising2D::getPosition(const uint particleIndex, uint &x, uint &y, uint &z) const
//{
//    x = 0;
//    y = 0;
//    z = 0;
//}
