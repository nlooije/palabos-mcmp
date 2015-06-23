
#include <cstdlib>
#include <iostream>
#include <math.h>

#include "palabos2D.h"
#include "palabos2D.hh"

// include modified libraries (set IncludePath in MakeFile)
#include "shanChenLattices2D.hm"
#include "isoThermalDynamics.hm"
#include "isoThermalDynamics.hhm"
#include "shanChenProcessor2D.hm"
#include "shanChenProcessor2D.hhm"
#include "interparticlePotential.hm"
#include "interparticlePotential.hhm"

using namespace plb;
using namespace std;

typedef double T;

#define DESCRIPTOR descriptors::GeneralForcedShanChenD2Q9Descriptor
#define DYNAMICS GeneralExternalMomentBGKdynamics<T, DESCRIPTOR> ( omega )
#define POTENTIAL interparticlePotential::PsiYuanSchaefer<T> ( Tr, eos )
#define PROCESSOR ShanChenSingleComponentProcessor2D_mod<T,DESCRIPTOR> ( G, new POTENTIAL, accel, rho_avg, kappa )
#define CASE RandomInitializer<T,DESCRIPTOR> ( rhoc, deltaRho )
// #define CASE DropletInitializer<T,DESCRIPTOR> ( rho_vap, rho_liq, iX0, iY0, radius, width )

// Interface for spinodal decomposition
template<typename T, template<typename U> class Descriptor>
class RandomInitializer : public OneCellIndexedFunctional2D<T,Descriptor> {
public:
    RandomInitializer(T rhoc_, T maxRho_) : rhoc(rhoc_), maxRho(maxRho_)
    { }
    RandomInitializer<T,Descriptor>* clone() const {
        return new RandomInitializer<T,Descriptor>(*this);
    }
    virtual void execute(plint iX, plint iY, Cell<T,Descriptor>& cell) const {
        T rand = (T)2*(T)random()/(T)RAND_MAX-(T)1;
        T rho = rhoc + rand * maxRho;
        Array<T,2> zeroVelocity (0.,0.);
        iniCellAtEquilibrium(cell, rho, zeroVelocity);
    }
private:
    T rhoc, maxRho;
};

// Interface for droplet initialization
template<typename T, template<typename U> class Descriptor>
class DropletInitializer : public OneCellIndexedFunctional2D<T,Descriptor> {
public:
    DropletInitializer(T rho_vap_, T rho_liq_, T iX0_, T iY0_, T radius_, T width_)
        : rho_vap(rho_vap_), rho_liq(rho_liq_), 
          iX0(iX0_), iY0(iY0_), radius(radius_), width(width_)
    { }
    DropletInitializer<T,Descriptor>* clone() const {
        return new DropletInitializer<T,Descriptor>(*this);
    }
    virtual void execute(plint iX, plint iY, Cell<T,Descriptor>& cell) const {
        T alpha = (T)0.5*(rho_vap + rho_liq); 
        T beta = (T)0.5*(rho_vap - rho_liq);
        T r = sqrt(((T)iX-iX0)*((T)iX-iX0)+((T)iY-iY0)*((T)iY-iY0));
        T rho = alpha + beta*tanh(-(r-radius)/width);
        Array<T,2> zeroVelocity (0.,0.);
        iniCellAtEquilibrium(cell, rho, zeroVelocity);
    }
private:
    T rho_vap, rho_liq;
    T iX0, iY0, radius, width;
};

template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice, pluint iter, pluint digits)
{
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, digits), 1);
    vtkOut.writeData<float>(*computeDensity(lattice), "density", 1);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", 1);
}

template<class BlockLatticeT>
void writeGif(BlockLatticeT& lattice, pluint iT, pluint digits)
{
    const plint imSize = 400;
    ImageWriter<T> imageWriter("leeloo.map");
    imageWriter.writeScaledGif(createFileName("rho", iT, digits),
                            *computeDensity(lattice),
                            imSize, imSize);
//    imageWriter.writeScaledGif(createFileName("umag", iT, digits),
//                            *computeVelocityNorm(lattice),
//                            imSize, imSize);
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    const pluint nx = 256+1;
    const pluint ny = nx;

    const T G = -1.0;
    const T Tr = 0.85;
    const std::string eos = "vdW"; // {vdW, RK, RKS, PR, CS}
    // const std::string forcing = "guo"; // {guo, shan-chen, edm, lycett-brown}
    const T kappa = 0;
    
    const pluint maxIter  = 1001;
    const pluint saveIter = 100;
    const pluint statIter = 100;

    // For spinodal decomposition:
    //     parameters depend on eos (see interparticlePotential.h)
    const T rhoc = 3.5; // critical density
    const T deltaRho = rhoc/(T)10; // rhoc/10

    // For droplet initialization:
    const T radius = (T)((ny-1)/4);
    const T width = (T)0.88;
    const T iX0 = (T)(1+(nx-1)/2);
    const T iY0 = (T)(1+(ny-1)/2);
    const T rho_vap = 1.17936;
    const T rho_liq = 6.38339;

    // van der Waals (Tr=0.85):
    // const T rhoc = 3.5; // critical density
    // const T rho_vap = 1.17936;
    // const T rho_liq = 6.38339;

    // Carnahan-Starling (Tr=0.85):
    // const T rhoc = 0.130444; // critical density
    // const T rho_vap = 0.0222236;
    // const T rho_liq = 0.27655;

    // requirements
    const T sigma = 0.1705; // see surface_tension/analysis.py
    const T diameter = (T)2.*radius;
    const T rho_ratio = (T)1.-rho_vap/rho_liq;

    // parameters
    const T Eo = 1.; // = gy*rho_liq*rho_ratio*D^2/sigma
    const T Mo = 1.; // = gy*rho_l^3*nu_l^4*rho_ratio/sigma^3
    const T Re = 1.; // = U*D/nu_l
    
    const T gy = sigma/(rho_liq*rho_ratio*pow(diameter,(T)2))*Eo;
    const T nu_liq = pow(pow(sigma,(T)3.)/(gy*pow(rho_liq,(T)3)*rho_ratio)*Mo,(T)0.25);
    const T omega = 1.;//(T)1./((T)0.5+(T)3.*nu_liq);
    const T umax = nu_liq/diameter*Re;
    
    pcout << "Eotvos: " << Eo << ", Morton: " << Mo << ", Reynolds: " << Re << endl; 
    pcout << "diameter: " << diameter << ", reducedT: " << Tr << ", eos: " << eos << endl;
    pcout << "rho_liq: " << rho_liq << ", rho_ratio: " << rho_ratio << ", sigma: " << sigma << endl;
    pcout << "gy: " << gy << ", tau: " << (T)1/omega << ", umax: " << umax << endl;  

    MultiBlockLattice2D<T, DESCRIPTOR> lattice ( nx, ny, new DYNAMICS );

    lattice.periodicity().toggleAll(true);

    applyIndexed(lattice, lattice.getBoundingBox(), new CASE );

    lattice.initialize();

    // Gravity force
    Array<T,2> accel(0.,0.);
    T rho_avg = computeAverage(*computeDensity(lattice));
    pcout << "rho_avg: " << rho_avg << endl;

    pcout << "Starting simulation" << endl;
    for (pluint iT=0; iT < maxIter; ++iT) {
        if (iT%statIter==0) {
            auto_ptr<MultiScalarField2D<T> > rho( computeDensity(lattice) );
            pcout << iT << ": " << endl;
            pcout << "Minimum density: " << computeMin(*rho) << endl;
            pcout << "Maximum density: " << computeMax(*rho) << endl;
            pcout << "Average denstiy: " << computeAverage(*rho)/rho_avg << endl;
        }
        if (iT%saveIter == 0) {
            pluint digits = 6;
            writeGif(lattice, iT, digits);
            // writeVTK(lattice, iT, digits);
        }
        applyProcessingFunctional( new PROCESSOR,
            lattice.getBoundingBox(), lattice );
        lattice.collideAndStream();
    }

    // write out density profile
    // plb_ofstream ofile(createFileName("./surface_tension/rho_profile_r=",radius,2).c_str());
    // ofile << setprecision(10) << *computeDensity(lattice, Box2D(0,nx-1,iY0,iY0)) << endl;
}

