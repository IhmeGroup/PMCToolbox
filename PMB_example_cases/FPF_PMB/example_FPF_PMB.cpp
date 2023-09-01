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
auto SiC_heatCond1 = std::make_shared<Cantera::Const1>(Cantera::Const1(1857.0));
auto SiC_heatCond2 = std::make_shared<Cantera::Pow1>(Cantera::Pow1(-0.5332));
auto YZA_heatCond1 = std::make_shared<Cantera::Const1>(Cantera::Const1(38.95));
auto YZA_heatCond2 = std::make_shared<Cantera::Pow1>(Cantera::Pow1(-0.35));

// function to create a list of three "FoamProperties" objects representing the three burner sections
std::vector<FoamProperties> createFoams()
{
    double emissivity       = 0.9;
    double transmissivity   = 0.6;

    double gapwidth = 0; // space between sections where all effective properties are linearly interpolated. Setting to e.g. 1e-4 m can sometimes improve convergence
/*
    FoamProperties YZA_40PPI;

    double YZA_40PPI_porosity           = 0.8249333333;
    double YZA_40PPI_specific_area      = 1591.815;

    YZA_40PPI.porosity                  = std::make_shared<Cantera::Const1>(Cantera::Const1(YZA_40PPI_porosity));
    YZA_40PPI.heat_conductivity.reset(new Cantera::Product1(*YZA_heatCond1, *YZA_heatCond2));
    YZA_40PPI.specific_area             = std::make_shared<Cantera::Const1>(Cantera::Const1(YZA_40PPI_specific_area));
    YZA_40PPI.hydraulic_diameter        = std::make_shared<Cantera::Const1>(Cantera::Const1(4.0*YZA_40PPI_porosity/YZA_40PPI_specific_area));
    YZA_40PPI.solidPhase                = Cantera::newSolution("data/YZA_new.yaml");
    YZA_40PPI.height                    = 2*inch - 0.5*gapwidth;
    YZA_40PPI.emissivity                = emissivity;
    YZA_40PPI.heatCond_eff_factor       = std::make_shared<Cantera::Const1>(Cantera::Const1(0.0716));
    YZA_40PPI.insulationTransmissivity  = transmissivity;
    YZA_40PPI.chemistryFactor           = 1;
    YZA_40PPI.tortuosity_factor         = std::make_shared<Cantera::Const1>(Cantera::Const1(1.343));
    YZA_40PPI.extinction_coefficient    = std::make_shared<Cantera::Const1>(Cantera::Const1(1340.4));

    FoamProperties SiC_3PPI;

    double SiC_3PPI_porosity            = 0.8623333333;
    double SiC_3PPI_specific_area       = 934.45;

    SiC_3PPI.porosity                   = std::make_shared<Cantera::Const1>(Cantera::Const1(SiC_3PPI_porosity));
    SiC_3PPI.heat_conductivity.reset(new Cantera::Product1(*SiC_heatCond1,*SiC_heatCond2));
    SiC_3PPI.specific_area              = std::make_shared<Cantera::Const1>(Cantera::Const1(SiC_3PPI_specific_area));
    SiC_3PPI.hydraulic_diameter         = std::make_shared<Cantera::Const1>(Cantera::Const1(4.0*SiC_3PPI_porosity/SiC_3PPI_specific_area));
    SiC_3PPI.solidPhase                 = Cantera::newSolution("data/SiC_new.yaml");
    SiC_3PPI.height                     = inch - gapwidth;
    SiC_3PPI.emissivity                 = emissivity;
    SiC_3PPI.heatCond_eff_factor        = std::make_shared<Cantera::Const1>(Cantera::Const1(0.0506));
    SiC_3PPI.insulationTransmissivity   = transmissivity;
    SiC_3PPI.chemistryFactor            = 1;
    SiC_3PPI.tortuosity_factor          = std::make_shared<Cantera::Const1>(Cantera::Const1(1.168));
    SiC_3PPI.extinction_coefficient     = std::make_shared<Cantera::Const1>(Cantera::Const1(525.95));
*/

    FoamProperties SiC_10PPI;

    double SiC_10PPI_specific_area      = 986.375;
    double SiC_10PPI_porosity           = 0.86;

    SiC_10PPI.porosity                  = std::make_shared<Cantera::Const1>(Cantera::Const1(SiC_10PPI_porosity));
    SiC_10PPI.heat_conductivity.reset(new Cantera::Product1(*SiC_heatCond1,*SiC_heatCond2));
    SiC_10PPI.specific_area             = std::make_shared<Cantera::Const1>(Cantera::Const1(SiC_10PPI_specific_area));
    SiC_10PPI.hydraulic_diameter        = std::make_shared<Cantera::Const1>(Cantera::Const1(4.0*SiC_10PPI_porosity/SiC_10PPI_specific_area));
    SiC_10PPI.solidPhase                = Cantera::newSolution("data/SiC_new.yaml");
    SiC_10PPI.height                    = 4*inch - 0.5*gapwidth;
    SiC_10PPI.emissivity                = emissivity;
    SiC_10PPI.heatCond_eff_factor       = std::make_shared<Cantera::Const1>(Cantera::Const1(0.0432));
    SiC_10PPI.insulationTransmissivity  = transmissivity;
    SiC_10PPI.chemistryFactor           = 1;
    SiC_10PPI.tortuosity_factor         = std::make_shared<Cantera::Const1>(Cantera::Const1(1.15014));
    SiC_10PPI.extinction_coefficient    = std::make_shared<Cantera::Const1>(Cantera::Const1(683.));

/*
    return  {
                YZA_40PPI,
                FoamProperties::createGap(gapwidth),
                SiC_3PPI,
                FoamProperties::createGap(gapwidth),
                SiC_10PPI
            };
*/
    return {SiC_10PPI};
}


int main ()
{
    auto sol = Cantera::newSolution("data/Arunthanayothin.xml");

    double diameter = 2*inch;
    FoamStack foam_stack(createFoams(), diameter);


    double T0 = 300.0;
    double Tamb = T0;
    double p0 = 101325.;

    double initialFlameProfileThickness = 0.01;
    std::size_t nInitialPoints = 10;
    double totalLength = foam_stack.totalHeight;
    std::string fuel = "NH3:0.85,O2:0.15";
    std::string oxidizer = "O2:0.21,N2:0.79";

    std::vector<double> phis{1.0};
    Cantera::vector_fp mdots{0.40};

    bool freeflame = true; // decide if this should be a matrix-stabilized flame or freely propagating flame
    bool radiation = true; // decide if you want to include radial and axial radiation in the simulation
    int iterations = 50;
    solveFlames(*sol, foam_stack, mdots, phis, fuel, oxidizer, T0, p0, Tamb,
        initialFlameProfileThickness, totalLength, nInitialPoints, 0.5*totalLength,
        freeflame,radiation,iterations,false,1);

    std::exit(EXIT_SUCCESS);
}

