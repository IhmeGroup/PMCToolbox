// ========================================================================
//  Thorsten Zirwes, Karlsruhe Institute of Technology
//  Guillaume Vignat, Stanford University
//  July 2023
//
//  If you use this code, please cite as:
//
//  T. Zirwes, G. Vignat, E.R. Toro, E. Boigne, K. Younes, D. Trimis, M. Ihme
//  Improving volume-averaged simulations of matrix-stabilized combustion through direct X-ray
//  muCT characterization: Application to NH3/H2-air combustion
//  Combustion and Flame, 2023
// ========================================================================

#include "example_PMB.h"

using namespace Cantera;

/*
 * Step 1: define the effective properties
 * Here, we assume three different sections comprising the burner
*/

double inch = 0.0254; // conversion factor between inches and meters


// heat conductivites of SiC and YZA
auto YZA_heatCond1 = std::make_shared<Cantera::Const1>(Cantera::Const1(38.95));
auto YZA_heatCond2 = std::make_shared<Cantera::Pow1>(Cantera::Pow1(-0.35));
double YZA_40PPI_porosity           = 0.8;

// function to create a list of three "FoamProperties" objects representing the three burner sections
std::vector<FoamProperties> createYZA()
{
    double emissivity       = 0.9;
    double transmissivity   = 0.6;

    double gapwidth = 0; // space between sections where all effective properties are linearly interpolated. Setting to e.g. 1e-4 m can sometimes improve convergence

    FoamProperties YZA_40PPI;

    double YZA_40PPI_specific_area      = 1591.815;

    YZA_40PPI.porosity                  = std::make_shared<Cantera::Const1>(Cantera::Const1(YZA_40PPI_porosity));
    YZA_40PPI.heat_conductivity.reset(new Cantera::Product1(*YZA_heatCond1, *YZA_heatCond2));
    YZA_40PPI.specific_area             = std::make_shared<Cantera::Const1>(Cantera::Const1(YZA_40PPI_specific_area));
    YZA_40PPI.hydraulic_diameter        = std::make_shared<Cantera::Const1>(Cantera::Const1(4.0*YZA_40PPI_porosity/YZA_40PPI_specific_area));
    YZA_40PPI.solidPhase                = Cantera::newSolution("data/YZA_new.yaml");
    YZA_40PPI.height                    = 0.1;
    YZA_40PPI.emissivity                = emissivity;
    YZA_40PPI.heatCond_eff_factor       = std::make_shared<Cantera::Const1>(Cantera::Const1(0.0716));
    YZA_40PPI.insulationTransmissivity  = transmissivity;
    YZA_40PPI.chemistryFactor           = 1;
    YZA_40PPI.tortuosity_factor         = std::make_shared<Cantera::Const1>(Cantera::Const1(1));
    YZA_40PPI.extinction_coefficient    = std::make_shared<Cantera::Const1>(Cantera::Const1(1340.4));

    return {YZA_40PPI};
}

int main()
{
    auto sol = Cantera::newSolution("data/Arunthanayothin.xml");

    double diameter = 2*inch;
    FoamStack foam_stack(createYZA(), diameter);


    double T0 = 300.0;
    double Tamb = T0;
    double p0 = 101325.;

    double initialFlameProfileThickness = 0.025;
    std::size_t nInitialPoints = 3;
    double totalLength = foam_stack.totalHeight;
    std::string fuel = "NH3";
    std::string oxidizer = "O2:0.21,N2:0.79";

    std::vector<double> phis{1.0};
    Cantera::vector_fp mdots{0.03*YZA_40PPI_porosity}; // initial guess for the FPF inlet mass flow

    bool freeflame = false; // decide if this should be a matrix-stabilized flame or freely propagating flame
    bool radiation = false; // decide if you want to include radial and axial radiation in the simulation
    bool deactivate_solid = true;
    int iterations = 20;

    solveFlames(*sol, foam_stack, mdots, phis, fuel, oxidizer, T0, p0, Tamb,
        initialFlameProfileThickness, totalLength, nInitialPoints, 1.1*initialFlameProfileThickness,
        freeflame,radiation,iterations,deactivate_solid,1);

    std::exit(EXIT_SUCCESS);
}
