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

#include "../example_PMB.h"

using namespace Cantera;

/*
 * Step 1: define the effective properties
 * Here, we assume three different sections comprising the burner
*/

double inch = 0.0254; // conversion factor between inches and meters

auto YZA_heatCond1 = std::make_shared<Cantera::Const1>(Cantera::Const1(38.95));
auto YZA_heatCond2 = std::make_shared<Cantera::Pow1>(Cantera::Pow1(-0.35));

// Linear increase of the pore diameter
auto pore_dia1 = std::make_shared<Cantera::Const1>(Cantera::Const1(0.00020));
auto pore_dia2 = std::make_shared<Cantera::Const1>(Cantera::Const1( (1.50-0.20) / 1000 / 4 / inch ));
auto pore_dia3 = std::make_shared<Cantera::Pow1>(Cantera::Pow1(1));
auto pore_dia_int = new Cantera::Product1(*pore_dia2, *pore_dia3);
// auto pore_dia_int = std::make_shared<Cantera::Product1>(Cantera::Product1(*pore_dia2, *pore_dia3));
auto pore_dia = new Cantera::Sum1(*pore_dia1, *pore_dia_int);
// auto pore_dia = std::make_shared<Cantera::Sum1>(Cantera::Sum1(*pore_dia1, *pore_dia_int));
auto Sv_den = std::make_shared<Cantera::Const1>(Cantera::Const1(4*0.8));
auto Sv = new Cantera::Ratio1(*Sv_den, *pore_dia);
auto extcoeff_den = std::make_shared<Cantera::Const1>(Cantera::Const1(3*(1-0.8)));
auto extcoeff = new Cantera::Ratio1(*extcoeff_den, *pore_dia);

double Length = 0.2;
double coeffs[] = {0.1,0.8/Length};
//double arr[] = {4,0.1};
auto porosity_function = Cantera::make_shared<Cantera::Poly1>(Cantera::Poly1(1, coeffs));

// function to create a list of three "FoamProperties" objects representing the three burner sections
// To add a region of heat exchange to a heat exchanger, add a fourth region in between the YZA and SiC 3 PPI
// in which we add some small heat losses
std::vector<FoamProperties> createFoams()
{
    double emissivity       = 0.0;
    double transmissivity   = 0.0;

    FoamProperties YZA;

    YZA.hydraulic_diameter        = std::make_shared<Cantera::Sum1>(*pore_dia);

    // make porosity change linearly from 0.1 to 0.9
    YZA.porosity                  = Cantera::make_shared<Cantera::Poly1>(Cantera::Poly1(1, coeffs)); // polynomial of degree 1
    YZA.heat_conductivity.reset(new Cantera::Product1(*YZA_heatCond1, *YZA_heatCond2));
    YZA.specific_area             = std::make_shared<Cantera::Ratio1>(*Sv);
    //YZA.volHeat_trans_coeff_HE    = std::make_shared<Cantera::Const1>(Cantera::Const1(0));
    YZA.solidPhase                = Cantera::newSolution("data/YZA_new.yaml");
    YZA.height                    = Length;
    YZA.emissivity                = emissivity;
    YZA.heatCond_eff_factor       = std::make_shared<Cantera::Const1>(Cantera::Const1(0.0716));
    YZA.insulationTransmissivity  = transmissivity;
    YZA.chemistryFactor           = 1;
    YZA.tortuosity_factor         = std::make_shared<Cantera::Const1>(Cantera::Const1(1.));
    YZA.extinction_coefficient    = std::make_shared<Cantera::Ratio1>(*extcoeff);

    return  {
                YZA
            };
}


int main ()
{
    auto sol = Cantera::newSolution("data/H2_K_10_12_0_OD.yaml");

    double diameter = 2*inch;
    FoamStack foam_stack(createFoams(), diameter);

    double T0 = 300.0;
    double Tamb = T0;
    double p0 = 101325.;

    double mdot_ref = 2.00000;

    double initialFlameProfileThickness = 0.005;
    double initialFlameLocation = 0.1;
    std::size_t nInitialPoints = 3;
    double totalLength = foam_stack.totalHeight;
    std::string fuel = "H2";
    std::string oxidizer = "O2:0.21,N2:0.79";

    Cantera::vector_fp phis{1.0};
    Cantera::vector_fp mdots{mdot_ref};

    bool freeflame = false; // decide if this should be a matrix-stabilized flame or freely propagating flame
    bool radiation = false; // decide if you want to include radial and axial radiation in the simulation
    int iterations = 300;

    solveFlames(*sol, foam_stack, mdots, phis, fuel, oxidizer, T0, p0, Tamb,
        initialFlameProfileThickness, totalLength, nInitialPoints, initialFlameLocation,
        freeflame,radiation,iterations,false,1);

    std::exit(EXIT_SUCCESS);
}
