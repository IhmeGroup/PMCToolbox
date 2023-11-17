// ========================================================================
//  Thorsten Zirwes, University of Stuttgart
//  Guillaume Vignat, Stanford University
//  November 2023
// ========================================================================


#include <cantera/ext/fmt/core.h>
#include <cantera/ext/fmt/printf.h>
#include <cantera/numerics/Func1.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include "cantera/base/Solution.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics.h"
#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/Boundary1D.h"
#include "cantera/oneD/Boundary1D_PMB.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/oneD/StFlow_PMB.h"
#include "cantera/transport.h"
#include <fstream>
#include <iomanip>
#include <type_traits>
#include "cantera/numerics/eigen_sparse.h"

using namespace Cantera;

const double T0 = 300.0;
const double Tsecond_in = T0;
const double Tamb = T0;
const double p0 = 101325.;
const double inch = 0.0254; // conversion factor between inches and meters

Eigen::SparseMatrix<double,Eigen::RowMajor> mat;
Eigen::VectorXd q;
Eigen::VectorXd radiation_rhs;
Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int> > solver;

template<typename T, typename U>
void solveTsolid(T& flow_first, T& flow_second, U& flame_first, U& flame_second,
        Cantera::vector_fp& prevMesh,
        Cantera::vector_fp& xc, Cantera::vector_fp& TsFixed)
{

    double m_Tsolid_inlet = T0;
    double m_T_amb = Tamb;
    size_t m_points = flow_first.nPoints() + flow_second.nPoints()-1;

    Cantera::vector_fp Tsolid_prev(m_points);
    Cantera::vector_fp m_effectiveSolidHeatConductivity(m_points);
    Cantera::vector_fp m_hv(m_points);
    Cantera::vector_fp m_scatteringAlbedo(m_points);
    Cantera::vector_fp m_extinctionCoefficient(m_points);
    Cantera::vector_fp m_emissivity(m_points);
    Cantera::vector_fp m_porosity(m_points);
    Cantera::vector_fp z(m_points);
    Cantera::vector_fp Tg(m_points);
    Cantera::vector_fp m_divq(m_points);
    Cantera::vector_fp m_rad_radial(m_points);

    for (size_t i=0; i!=flow_first.nPoints(); ++i)
    {
        m_effectiveSolidHeatConductivity[i] = flow_first.m_effectiveSolidHeatConductivity[i];
        m_hv[i] = flow_first.m_hv[i];
        m_scatteringAlbedo[i] = flow_first.m_scatteringAlbedo[i];
        m_extinctionCoefficient[i] = flow_first.m_extinctionCoefficient[i];
        //m_divq[i] = flow_first.m_divq[i];
        m_rad_radial[i] = flow_first.m_rad_radial[i];
        m_porosity[i] = flow_first.m_porosity[i];
        m_emissivity[i] = flow_first.m_emissivity[i];
        z[i] = flow_first.grid(i);
        Tg[i] = flame_first.value(1,flow_first.componentIndex("T"),i);
        Tsolid_prev[i] = flame_first.value(1,flow_first.componentIndex("Tsolid"),i);
    }
    for (size_t i=1; i!=flow_second.nPoints(); ++i)
    {
        m_effectiveSolidHeatConductivity[flow_first.nPoints() -1 + i] = flow_second.m_effectiveSolidHeatConductivity[i];
        m_hv[flow_first.nPoints() -1 + i] = flow_second.m_hv[i];
        m_scatteringAlbedo[flow_first.nPoints() -1 + i] = flow_second.m_scatteringAlbedo[i];
        m_extinctionCoefficient[flow_first.nPoints() -1 + i] = flow_second.m_extinctionCoefficient[i];
        //m_divq[flow_first.nPoints() -1 + i] = flow_second.m_divq[i];
        m_rad_radial[flow_first.nPoints() -1 + i] = flow_second.m_rad_radial[i];
        m_porosity[flow_first.nPoints() -1 + i] = flow_second.m_porosity[i];
        m_emissivity[flow_first.nPoints() -1 + i] = flow_second.m_emissivity[i];
        z[flow_first.nPoints() -1 + i] = flow_second.grid(i);
        Tg[flow_first.nPoints() -1 + i] = flame_second.value(1,flow_second.componentIndex("T"),i);
        Tsolid_prev[flow_first.nPoints() -1 + i] = flame_second.value(1,flow_second.componentIndex("Tsolid"),i);
    }

    {
        bool recalcMatrix = false;
        if (prevMesh.size() != m_points)
        {
            recalcMatrix = true;
        }
        if (!recalcMatrix)
        {
            for(std::size_t i=0; i!=m_points; ++i)
            {
                if (std::abs(z[i] - prevMesh[i]) > 1e-8)
                {
                    recalcMatrix = true;
                    break;
                }
            }
        }
        if (recalcMatrix)
        {
            prevMesh.resize(m_points);
            for(std::size_t i=0; i!=m_points; ++i)
                prevMesh[i] = z[i];
        }

        constexpr int non_zeros_per_col = 3;

        if (recalcMatrix)
        {
            mat = Eigen::SparseMatrix<double,Eigen::RowMajor>(2*m_points,2*m_points);
            mat.reserve(Eigen::Matrix<int, Eigen::Dynamic, 1>::Constant(2*m_points, non_zeros_per_col));

            mat.insert(0,0) = 1.0;
            for (int i=1; i!=static_cast<int>(m_points); ++i)
            {
                mat.insert(i, i-1) = -1.0;
                mat.insert(i, i) = 1.0+(z[i]-z[i-1])*m_extinctionCoefficient[i]*(2.0-m_scatteringAlbedo[i]);
                mat.insert(i, i+m_points) = -(z[i]-z[i-1])*m_extinctionCoefficient[i]*m_scatteringAlbedo[i];
            }
            for (int j=m_points; j!=2*static_cast<int>(m_points)-1; ++j)
            {
                int i = j-m_points;
                mat.insert(j, j-m_points) = -(z[i+1]-z[i])*m_extinctionCoefficient[i]*m_scatteringAlbedo[i];
                mat.insert(j, j) = 1.0+(z[i+1]-z[i])*m_extinctionCoefficient[i]*(2.0-m_scatteringAlbedo[i]);
                mat.insert(j, j+1) = -1.0;
            }
            mat.insert(2*m_points-1,2*m_points-1) = 1.0;
            mat.makeCompressed();

            q.resize(2*m_points);
            radiation_rhs.resize(2*m_points);
            radiation_rhs(0) = StefanBoltz*std::pow(m_T_amb,4.0);
            radiation_rhs(2*m_points-1) = radiation_rhs(0);

            solver.compute(mat);
            auto res = solver.info();
            if(res != Eigen::Success)
            {
                std::cout<<"ERROR: eigen stage 1 fail"<<std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (res == Eigen::NumericalIssue)
            {
                std::cout<<"ERROR: eigen stage 1 fail - numerics"<<std::endl;
                std::exit(EXIT_FAILURE);

            }
        }

        for (int i=1; i!=static_cast<int>(m_points); ++i)
            radiation_rhs(i) = (z[i]-z[i-1])*2.0*m_extinctionCoefficient[i]*(1.0-m_scatteringAlbedo[i])*StefanBoltz*std::pow(Tsolid_prev[i], 4.0);
        for (int j=m_points; j!=2*static_cast<int>(m_points)-1; ++j)
        {
            int i = j-m_points;
            radiation_rhs(j) = (z[i+1]-z[i])*2.0*m_extinctionCoefficient[i]*(1.0-m_scatteringAlbedo[i])*StefanBoltz*std::pow(Tsolid_prev[i], 4.0);
        }

        q = solver.solve(radiation_rhs);
        auto res = solver.info();
        if(res != Eigen::Success)
        {
            std::cout<<"ERROR: eigen stage 1 fail"<<std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (res == Eigen::NumericalIssue)
        {
            std::cout<<"ERROR: eigen stage 1 fail - numerics"<<std::endl;
            std::exit(EXIT_FAILURE);

        }

        for (int i=0; i!=static_cast<int>(m_points); ++i)
            m_divq[i] = 2*m_extinctionCoefficient[i]*(1.0-m_scatteringAlbedo[i])*(2.0*StefanBoltz*std::pow(Tsolid_prev[i], 4.0) - (q[i]+q[i+m_points]));
    }

    // =========================================================================================================================
    Eigen::SparseMatrix<double,Eigen::RowMajor> Tsolid_mat = Eigen::SparseMatrix<double,Eigen::RowMajor>(m_points,m_points);
    constexpr int non_zeros_per_col = 3;
    Eigen::VectorXd Tsolid_q;
    Eigen::VectorXd Tsolid_rhs;
    Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int> > Tsolid_solver;
    Tsolid_mat.reserve(Eigen::Matrix<int, Eigen::Dynamic, 1>::Constant(m_points, non_zeros_per_col));
    Tsolid_mat.insert(0,0) = 1.0;


    for (int i=1; i!=static_cast<int>(m_points-1); ++i)
    {

        double lambda_mid_left = 0.5*(m_effectiveSolidHeatConductivity[i-1]+m_effectiveSolidHeatConductivity[i]);
        double lambda_mid_right = 0.5*(m_effectiveSolidHeatConductivity[i+1]+m_effectiveSolidHeatConductivity[i]);

        double c1 = 2.0*lambda_mid_left /((z[i+1]-z[i-1])*(z[i]-z[i-1]));
        double c2 = 2.0*lambda_mid_right/((z[i+1]-z[i-1])*(z[i+1]-z[i]));

        Tsolid_mat.insert(i, i-1) = -c1;
        Tsolid_mat.insert(i, i) = m_hv[i] + c1 + c2;
        Tsolid_mat.insert(i, i+1) = -c2;
    }

    /*
    if (m_TsolidZeroGradient)
    {
        Tsolid_mat.insert(m_points-1,m_points-1) = 1;
        Tsolid_mat.insert(m_points-1,m_points-2) = -1;
    }
    else
    */
    {
        double lambda_end_mid = 0.5*(m_effectiveSolidHeatConductivity[m_points-1]+m_effectiveSolidHeatConductivity[m_points-2])/(z[m_points-1]-z[m_points-2]);
        Tsolid_mat.insert(m_points-1,m_points-1) = lambda_end_mid + m_emissivity[m_points-1]*Cantera::StefanBoltz*(1.-m_porosity[m_points-1])*Tsolid_prev[m_points-1]*Tsolid_prev[m_points-1]*Tsolid_prev[m_points-1];
        Tsolid_mat.insert(m_points-1,m_points-2) = -lambda_end_mid;
    }

    Tsolid_mat.makeCompressed();
    Tsolid_q.resize(m_points);
    Tsolid_rhs.resize(m_points);
    Tsolid_rhs(0) = m_Tsolid_inlet;

    Tsolid_solver.compute(Tsolid_mat);
    auto res = Tsolid_solver.info();
    if(res != Eigen::Success)
    {
        std::cout<<"ERROR: eigen stage 1 fail"<<std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (res == Eigen::NumericalIssue)
    {
        std::cout<<"ERROR: eigen stage 1 fail - numerics"<<std::endl;
        std::exit(EXIT_FAILURE);
    }

    // modify rhs

    for (int i=1; i!=static_cast<int>(m_points-1); ++i)
        Tsolid_rhs(i) = m_hv[i] * Tg[i] - m_divq[i] - m_rad_radial[i];

    /*
    if (m_TsolidZeroGradient)
    {
        Tsolid_rhs(m_points-1) = 0;//Tsolid_rhs(m_points-2);
    }
    else
    */
    {
        Tsolid_rhs(m_points-1) = m_emissivity[m_points-1]*StefanBoltz*(1.-m_porosity[m_points-1])*m_T_amb*m_T_amb*m_T_amb*m_T_amb;
    }

    Tsolid_q = Tsolid_solver.solve(Tsolid_rhs);
    auto res2 = Tsolid_solver.info();
    if(res2 != Eigen::Success)
    {
        std::cout<<"ERROR: eigen stage 2 Tsolid fail"<<std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (res2 == Eigen::NumericalIssue)
    {
        std::cout<<"ERROR: eigen stage 2 Tsolid fail - numerics"<<std::endl;
        std::exit(EXIT_FAILURE);
    }

    TsFixed.resize(m_points);
    xc.resize(m_points);

    double m_relax = 1.0;//0.7;
    for (int i=0; i!=static_cast<int>(m_points); ++i)
    {
        TsFixed[i] = m_relax*Tsolid_q(i) + (1.-m_relax)*Tsolid_prev[i];
        xc[i] = z[i];
    }
}



template<typename F>
double bisect(F& f, double a, double b)
{
    if (f(a) * f(b) >= 0)
    {
        std::cout<<"ERROR in bisect!"<<std::endl;;
        std::exit(EXIT_FAILURE);;
    }
    double c = a;
    while ((b-a) >= 1e-8)
    {
        c = 0.5*(a+b);
        if (std::abs(f(c)) < 1e-9)
            break;
        else if (f(c)*f(a) < 0.)
            b = c;
        else
            a = c;
    }
    return c;
}

// heat conductivites of SiC and YZA
auto SiC_heatCond1 = std::make_shared<Cantera::Const1>(Cantera::Const1(1857.0));
auto SiC_heatCond2 = std::make_shared<Cantera::Pow1>(Cantera::Pow1(-0.5332));
auto YZA_heatCond1 = std::make_shared<Cantera::Const1>(Cantera::Const1(38.95));
auto YZA_heatCond2 = std::make_shared<Cantera::Pow1>(Cantera::Pow1(-0.35));

double emissivity       = 0.9;
double transmissivity   = 0.6;

// function to create a list of three "FoamProperties" objects representing the three burner sections
std::vector<FoamProperties> createFoams()
{
    double gapwidth = 0; // space between sections where all effective properties are linearly interpolated. Setting to e.g. 1e-4 m can sometimes improve convergence

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


    FoamProperties SiC_10PPI;

    double SiC_10PPI_specific_area      = 986.375;
    double SiC_10PPI_porosity           = 0.86;

    SiC_10PPI.porosity                  = std::make_shared<Cantera::Const1>(Cantera::Const1(SiC_10PPI_porosity));
    SiC_10PPI.heat_conductivity.reset(new Cantera::Product1(*SiC_heatCond1,*SiC_heatCond2));
    SiC_10PPI.specific_area             = std::make_shared<Cantera::Const1>(Cantera::Const1(SiC_10PPI_specific_area));
    SiC_10PPI.hydraulic_diameter        = std::make_shared<Cantera::Const1>(Cantera::Const1(4.0*SiC_10PPI_porosity/SiC_10PPI_specific_area));
    SiC_10PPI.solidPhase                = Cantera::newSolution("data/SiC_new.yaml");
    SiC_10PPI.height                    = inch - 0.5*gapwidth;
    SiC_10PPI.emissivity                = emissivity;
    SiC_10PPI.heatCond_eff_factor       = std::make_shared<Cantera::Const1>(Cantera::Const1(0.0432));
    SiC_10PPI.insulationTransmissivity  = transmissivity;
    SiC_10PPI.chemistryFactor           = 1;
    SiC_10PPI.tortuosity_factor         = std::make_shared<Cantera::Const1>(Cantera::Const1(1.15014));
    SiC_10PPI.extinction_coefficient    = std::make_shared<Cantera::Const1>(Cantera::Const1(683.));

    return  {
                YZA_40PPI,
                FoamProperties::createGap(gapwidth),
                SiC_3PPI,
                FoamProperties::createGap(gapwidth),
                SiC_10PPI
            };
}
std::vector<FoamProperties> createFoamsSecondStage()
{
    double gapwidth = 0; // space between sections where all effective properties are linearly interpolated. Setting to e.g. 1e-4 m can sometimes improve convergence


    FoamProperties SiC_3PPI;

    double SiC_3PPI_porosity            = 0.8623333333;
    double SiC_3PPI_specific_area       = 934.45;

    SiC_3PPI.porosity                   = std::make_shared<Cantera::Const1>(Cantera::Const1(SiC_3PPI_porosity));
    SiC_3PPI.heat_conductivity.reset(new Cantera::Product1(*SiC_heatCond1,*SiC_heatCond2));
    SiC_3PPI.specific_area              = std::make_shared<Cantera::Const1>(Cantera::Const1(SiC_3PPI_specific_area));
    SiC_3PPI.hydraulic_diameter         = std::make_shared<Cantera::Const1>(Cantera::Const1(4.0*SiC_3PPI_porosity/SiC_3PPI_specific_area));
    SiC_3PPI.solidPhase                 = Cantera::newSolution("data/SiC_new.yaml");
    SiC_3PPI.height                     = inch - 0.5*gapwidth;
    SiC_3PPI.emissivity                 = emissivity;
    SiC_3PPI.heatCond_eff_factor        = std::make_shared<Cantera::Const1>(Cantera::Const1(0.0506));
    SiC_3PPI.insulationTransmissivity   = transmissivity;
    SiC_3PPI.chemistryFactor            = 1;
    SiC_3PPI.tortuosity_factor          = std::make_shared<Cantera::Const1>(Cantera::Const1(1.168));
    SiC_3PPI.extinction_coefficient     = std::make_shared<Cantera::Const1>(Cantera::Const1(525.95));


    FoamProperties SiC_10PPI;

    double SiC_10PPI_specific_area      = 986.375;
    double SiC_10PPI_porosity           = 0.86;

    SiC_10PPI.porosity                  = std::make_shared<Cantera::Const1>(Cantera::Const1(SiC_10PPI_porosity));
    SiC_10PPI.heat_conductivity.reset(new Cantera::Product1(*SiC_heatCond1,*SiC_heatCond2));
    SiC_10PPI.specific_area             = std::make_shared<Cantera::Const1>(Cantera::Const1(SiC_10PPI_specific_area));
    SiC_10PPI.hydraulic_diameter        = std::make_shared<Cantera::Const1>(Cantera::Const1(4.0*SiC_10PPI_porosity/SiC_10PPI_specific_area));
    SiC_10PPI.solidPhase                = Cantera::newSolution("data/SiC_new.yaml");
    SiC_10PPI.height                    = inch - gapwidth;
    SiC_10PPI.emissivity                = emissivity;
    SiC_10PPI.heatCond_eff_factor       = std::make_shared<Cantera::Const1>(Cantera::Const1(0.0432));
    SiC_10PPI.insulationTransmissivity  = transmissivity;
    SiC_10PPI.chemistryFactor           = 1;
    SiC_10PPI.tortuosity_factor         = std::make_shared<Cantera::Const1>(Cantera::Const1(1.15014));
    SiC_10PPI.extinction_coefficient    = std::make_shared<Cantera::Const1>(Cantera::Const1(683.));

    FoamProperties SiC_20PPI;

    double SiC_20PPI_porosity            = 0.86;
    double SiC_20PPI_specific_area       = 1480.;

    SiC_20PPI.porosity                   = std::make_shared<Cantera::Const1>(Cantera::Const1(SiC_20PPI_porosity));
    SiC_20PPI.heat_conductivity.reset(new Cantera::Product1(*SiC_heatCond1,*SiC_heatCond2));
    SiC_20PPI.specific_area              = std::make_shared<Cantera::Const1>(Cantera::Const1(SiC_20PPI_specific_area));
    SiC_20PPI.hydraulic_diameter         = std::make_shared<Cantera::Const1>(Cantera::Const1(4.0*SiC_20PPI_porosity/SiC_20PPI_specific_area));
    SiC_20PPI.solidPhase                 = Cantera::newSolution("data/SiC_new.yaml");
    SiC_20PPI.height                     = inch - 0.5*gapwidth;
    SiC_20PPI.emissivity                 = emissivity;
    SiC_20PPI.heatCond_eff_factor        = std::make_shared<Cantera::Const1>(Cantera::Const1(1.0-0.86));
    SiC_20PPI.insulationTransmissivity   = transmissivity;
    SiC_20PPI.chemistryFactor            = 1;
    SiC_20PPI.tortuosity_factor          = std::make_shared<Cantera::Const1>(Cantera::Const1(1.0));
    SiC_20PPI.extinction_coefficient     = std::make_shared<Cantera::Const1>(Cantera::Const1(860));


    return  {
                SiC_3PPI,
                FoamProperties::createGap(gapwidth),
                SiC_10PPI,
                FoamProperties::createGap(gapwidth),
                SiC_20PPI
            };
}


int main ()
{
    using namespace Cantera;

    // some general settings
    auto solu = Cantera::newSolution("data/Arunthanayothin.xml");
    auto& sol = *solu;
    double diameter = 2*inch;
    FoamStack foam_stack(createFoams(), diameter);
    FoamStack foam_stack_secondStage(createFoamsSecondStage(), diameter);
    // settings for the first stage
    double initialFlameProfileThickness = 0.025;
    std::size_t nInitialPoints = 3;
    double totalLength = foam_stack.totalHeight;
    double totalLength_second = foam_stack_secondStage.totalHeight;
    double expected_flame_location = 0.5*totalLength;
    bool deactivateSolid = false;
    int loglevel = 1;
    Cantera::vector_fp initial_Ts_x = {0., 0.0498, 0.0518, 1.5}; // locations for initial temperature values
    Cantera::vector_fp initial_Ts   = {300., 300., 1600., 1200.}; // intital solid temperatures

    // ================================================
    // ==================  INPUT ======================

    // conditions first stage
    std::string fuel = "NH3";
    std::string oxidizer = "O2:0.21,N2:0.79";
    std::vector<double> phis{1.3};
    Cantera::vector_fp mdots{0.30};

    //conditions second stage:
    double phi_p = phis[0];
    double mdot_p = mdots[0];
    double phi_g = 0.95;
    double mdot_g = -1.;
    double phi_s = 0.5; // asummes hydrogen/air

    // =============================================



    bool second_stage_pure_air = false;
    if (phi_s < 0.01)
        second_stage_pure_air = true;

    std::cout<<"HERE1"<<std::endl;

    if (second_stage_pure_air)
    {
        phi_s = 0.;
    }
    // do some input sanity checks
    if (second_stage_pure_air && phi_s > 0)
    {
        std::cout<<"ERROR! Only one of {phi_s, second_stage_pure_air} can be set!"<<std::endl;
        std::exit(EXIT_FAILURE);
    }
    if ( (phi_g > 0) + (mdot_g>0) /*+ (mdot_p>0)*/ != 1)
    {
        std::cout<<"ERROR! Only one of {phi_g,mdot_g} can be set!"<<std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::cout<<"HERE2"<<std::endl;
    auto solu_del = Cantera::newSolution("data/Arunthanayothin.xml");
    auto& sol_del = *solu_del;
    auto& gas_del = *sol.thermo();
    gas_del.setState_TP(T0, p0);
    gas_del.setEquivalenceRatio(phi_p,fuel,oxidizer);
    Cantera::vector_fp Y_first(gas_del.nSpecies());
    gas_del.getMassFractions(Y_first.data());
    if (second_stage_pure_air)
        gas_del.setState_TPX(T0,p0,"O2:0.21,N2:0.79");
    else
        gas_del.setEquivalenceRatio(phi_s,"H2:1","O2:0.21,N2:0.79");
    vector_fp Y_second(gas_del.nSpecies());
    gas_del.getMassFractions(Y_second.data());
    gas_del.setState_TPY(Tsecond_in, p0, Y_second.data());
    double hsecond_in = gas_del.enthalpy_mass();

    std::cout<<"HERE3"<<std::endl;
    if (phi_g < 0.)
    {
        Cantera::vector_fp Y_total(gas_del.nSpecies());
        for(size_t i=0; i!=gas_del.nSpecies(); ++i)
        {
            Y_total[i] = mdot_p * Y_first[i] + (mdot_g-mdot_p)*Y_second[i];
        }
        gas_del.setState_TPY(T0,p0,Y_total.data());
        phi_g = gas_del.equivalenceRatio();
    }

    if (mdot_g < 0)
    {
        auto get_phi_g = [&](double mdot_gg)
        {
            Cantera::vector_fp Y_total(gas_del.nSpecies());
            for(size_t i=0; i!=gas_del.nSpecies(); ++i)
                Y_total[i] = mdot_p * Y_first[i] + (mdot_gg-mdot_p)*Y_second[i];

            gas_del.setState_TPY(T0,p0,Y_total.data());
            return gas_del.equivalenceRatio() - phi_g;
        };
        mdot_g = bisect(get_phi_g, mdot_p, 10);
    }

    if (mdot_p < 0.0)
    {
        auto get_phi_g = [&](double mdot_p)
        {
            Cantera::vector_fp Y_total(gas_del.nSpecies());
            for(size_t i=0; i!=gas_del.nSpecies(); ++i)
                Y_total[i] = mdot_p * Y_first[i] + (mdot_g-mdot_p)*Y_second[i];

            gas_del.setState_TPY(T0,p0,Y_total.data());
            return gas_del.equivalenceRatio() - phi_g;
        };
        mdot_g = bisect(get_phi_g, 0, mdot_g);
    }

    double mdot_s = mdot_g - mdot_p;
    std::cout<<"Second stage info:\n";
    std::cout<<"mdot_g: "<<mdot_g<<'\n';
    std::cout<<"mdot_s: "<<mdot_s<<'\n';
    std::cout<<"phi_s: "<<phi_s<<'\n';
    std::cout<<"phi_g: "<<phi_g<<'\n';

    bool freeflame = false; // decide if this should be a matrix-stabilized flame or freely propagating flame
    bool radiation = true; // decide if you want to include radial and axial radiation in the simulation
    int iterations = 500; // maximum number of iterations between gas and solid phase
    double residual = 0.1; // stop iterating once the residuals of both gas and solid temperature are below this value
    double residual_second = 0.01;
    double final_residual = 0.01;
    double relax = 1.0; // explicit under-relaxation factor

    auto& gas = *sol.thermo();
    size_t nsp = gas.nSpecies();

    // setup the first stage
    double pressure = p0;
    double mdot0 = mdots[0];
    double phi = phis[0];
    gas.setEquivalenceRatio(phi, fuel, oxidizer);
    gas.setState_TP(T0, pressure);
    Cantera::vector_fp x(nsp);
    gas.getMoleFractions(x.data());
    double rho_in = gas.density();
    Cantera::vector_fp yin(nsp);
    gas.getMassFractions(&yin[0]);
    gas.equilibrate("HP");
    double Tb = gas.temperature();
    Cantera::vector_fp yout(gas.nSpecies());
    gas.getMassFractions(&yout[0]);
    gas.setState_TP(Tb, pressure);
    gas.getMassFractions(&yout[0]);
    double rho_out = gas.density();
    double Tad = gas.temperature();
    gas.setEquivalenceRatio(phi, fuel, oxidizer);
    gas.setState_TP(T0, pressure);
    Cantera::PorousFlow flow(foam_stack, &gas);
    Cantera::vector_fp z(nInitialPoints+2);
    z[0]=0;
    auto dz = initialFlameProfileThickness / (nInitialPoints - 1);
    for (std::size_t iz = 0; iz < nInitialPoints; iz++)
        z[iz + 1] = iz*dz + expected_flame_location - initialFlameProfileThickness*0.5;
    z[nInitialPoints+1] = totalLength;
    flow.setupGrid(z.size(), &z[0]);
    std::unique_ptr<Cantera::Transport> trmix(Cantera::newTransportMgr("Mix", sol.thermo().get()));
    flow.setTransport(*trmix);
    flow.setKinetics(*sol.kinetics());
    flow.setPressure(pressure);
    flow.set_Tamb(Tamb);
    flow.set_Tsolid_inlet(T0);
    flow.setSolid(!deactivateSolid);
    flow.setSolidAxialRadiation(radiation);
    flow.setSolidRadialRadiation(radiation);
    flow.setAxisymmetricFlow();
    Cantera::InletWithSolid1D inlet;
    inlet.setFoams(&flow.foams);
    inlet.setMoleFractions(x.data());
    double uin = mdot0/(rho_in*foam_stack.get_local_foam_properties(0.0, T0).porosity);
    inlet.setMdot(mdot0);
    inlet.setTemperature(T0);
    inlet.setTsolid(T0);
    Cantera::OutletWithSolid1D outlet;
    outlet.setFoams(&flow.foams);
    outlet.setTamb(Tamb);
    std::vector<Cantera::Domain1D*> domains { &inlet, &flow, &outlet };
    Cantera::Sim1D flame(domains);
    flame.setGridMin(1,1e-6);
    Cantera::vector_fp locs{0.0, (expected_flame_location-0.5*initialFlameProfileThickness)/totalLength,
        (expected_flame_location+0.5*initialFlameProfileThickness)/totalLength, 1.0};
    double porosity_out = foam_stack.get_local_foam_properties(totalLength, Tad).porosity;
    double uout = mdot0/(rho_out*porosity_out);
    Cantera::vector_fp value{uin, uin, uout, uout};
    flame.setInitialGuess("velocity",locs,value);
    value = {T0, T0, Tad, Tad};
    flame.setInitialGuess("T",locs,value);
    for (size_t i=0; i<nsp; i++)
    {
        value = {yin[i], yin[i], yout[i], yout[i]};
        flame.setInitialGuess(gas.speciesName(i),locs,value);
    }
    value = {T0, T0, Tad, Tad};
    flame.setInitialGuess("Tsolid",locs,value);
    inlet.setMoleFractions(x.data());
    inlet.setMdot(mdot0);
    inlet.setTemperature(T0);
    inlet.setTsolid(T0);
    int flowdomain = 1;
    flow.set_initial_guess_Ts(initial_Ts_x, initial_Ts);
    flow.set_relaxation(relax);
    flow.set_fixedTemperature(true);
    flow.solveEnergyEqn();

    auto setTx = [&flow,&flame,flowdomain](Cantera::vector_fp& x, Cantera::vector_fp& T, Cantera::vector_fp& Ts)
    {
        T.resize(flow.nPoints());
        Ts.resize(flow.nPoints());
        x.resize(flow.nPoints());
        for (std::size_t i=0; i!=flow.nPoints(); ++i)
        {
            x[i] = flow.grid(i);
            Ts[i] = flame.value(flowdomain,flow.componentIndex("Tsolid"),i);
            T[i] = flame.value(flowdomain,flow.componentIndex("T"),i);
        }
    };

    auto interpolate_T = [](const Cantera::vector_fp& vxnew, const Cantera::vector_fp& vxold, const Cantera::vector_fp& Told, Cantera::vector_fp& Tnew)
    {
        Tnew.resize(vxnew.size());
        size_t pos = 0;
        for(std::size_t p=0; p!=vxnew.size(); ++p)
        {
            double xnew = vxnew[p];
            while(true)
            {
                double xold = vxold[pos];
                if (std::abs(xold-xnew)<1e-12)
                {
                    Tnew[p] = Told[pos];
                    break;
                }
                else if (xold < xnew && pos == vxold.size()-1)
                {
                    Tnew[p] = Told[pos];
                    break;
                }
                else if (xold <= xnew && vxold[pos+1] >= xnew)
                {
                    Tnew[p] = Told[pos] + (xnew-xold)*(Told[pos+1]-Told[pos])/(vxold[pos+1]-vxold[pos]);
                    break;
                }
                else
                    ++pos;
            }
        }
    };

    auto L2 = [](const Cantera::vector_fp& a, const Cantera::vector_fp& b)
    {
        double sum = 0.;
        for (size_t i=0; i!=a.size(); ++i)
            sum += (a[i]-b[i]) * (a[i]-b[i]);
        return std::sqrt(sum/a.size());
    };

    // solve the first stage
    {
        flow.setTransientTolerances(3e-2, 1e-6);
        flow.setSteadyTolerances(3e-2, 1e-6);
        flame.setRefineCriteria(flowdomain, 3, 0.5, 0.5, -0.05);
        flame.solve(loglevel,true);
        double ratio = 2.0;
        double slope = 0.15;
        double curve = 0.25;
        double prune = 0.01;
        flame.setRefineCriteria(flowdomain, ratio, slope, curve, prune);
        Cantera::vector_fp Xp(gas.nSpecies());
        gas.setEquivalenceRatio(phi, fuel, oxidizer);
        gas.setState_TP(T0, pressure);
        Cantera::vector_fp x(nsp);
        gas.getMoleFractions(x.data());
        inlet.setMoleFractions(x.data());
        inlet.setMdot(mdot0);
        flame.solve(loglevel,true);
        flow.setTransientTolerances(1e-4, 1e-9);
        flow.setSteadyTolerances(1e-4, 1e-9);

        flame.solve(loglevel,true);

        flow.setTransientTolerances(1e-7, 1e-13);
        flow.setSteadyTolerances(1e-7, 1e-13);

        Cantera::vector_fp gasT_prev;
        Cantera::vector_fp gasTs_prev;
        Cantera::vector_fp gasx_prev;
        Cantera::vector_fp gasT_curr;
        Cantera::vector_fp gasTs_curr;
        Cantera::vector_fp gasx_curr;
        Cantera::vector_fp T_tmp;

        setTx(gasx_prev, gasT_prev, gasTs_prev);
        for (int i=0; i!=iterations; ++i)
        {
            flow.set_firstSolve(true);
            flame.solve(loglevel,true);

            setTx(gasx_curr, gasT_curr, gasTs_curr);
            std::cout<<"Coupling gas and solid: "<<i+1<<"/"<<iterations<<'\n';

            // interpolate old results onto grid of new solution
            interpolate_T(gasx_curr, gasx_prev, gasT_prev, T_tmp);
            double L2gas = L2(gasT_curr, T_tmp);
            std::cout<<"initial solution stage 1: L2 Tg: "<<L2gas<<'\n';
            interpolate_T(gasx_curr, gasx_prev, gasTs_prev, T_tmp);
            double L2solid = L2(gasTs_curr, T_tmp);
            std::cout<<"initial solution stage 1: L2 Ts: "<<L2solid<<'\n';

            if (L2gas < residual && L2solid < residual)
            {
                std::cout<<'\n'<<"Phases converged!"<<'\n';
                break;
            }
            gasx_curr.swap(gasx_prev);
            gasT_curr.swap(gasT_prev);
            gasTs_curr.swap(gasTs_prev);
        }

        std::cout<<"success"<<std::endl;
    }

    //=====================================================================
    // compute the gas-phase inlet conditions for the second stage

    Cantera::vector_fp Yfirst_out(nsp);
    for(size_t k=0; k!=nsp; ++k)
        Yfirst_out[k] = flame.value(flowdomain,flow.componentIndex(gas.speciesName(k)),flow.nPoints()-1);
    double Tfirst_out = flame.value(flowdomain,flow.componentIndex("T"),flow.nPoints()-1);
    gas_del.setState_TPY(Tfirst_out,p0,Yfirst_out.data());
    double hfirst_out = gas_del.enthalpy_mass();
    double Tsfirst_out = flame.value(flowdomain,flow.componentIndex("Tsolid"),flow.nPoints()-1);

    Cantera::vector_fp YSecond_in(nsp);
    for(size_t k=0; k!=nsp; ++k)
        YSecond_in[k] = mdot_p * Yfirst_out[k] + mdot_s * Y_second[k];
    gas_del.setState_TPY(300.,p0,YSecond_in.data());
    gas_del.getMassFractions(YSecond_in.data());
    double hsecond_in_combined = mdot_p/mdot_g * hfirst_out + mdot_s/mdot_g*hsecond_in;
    gas_del.setState_HP(hsecond_in_combined, p0);
    double Tsecond_in_combined = gas_del.temperature();
    double rhosecond_in_combined = gas_del.density();
    std::cout<<"Inlet conditions for second stage:"<<std::endl;
    std::cout<<"Tfirst out: "<<Tfirst_out<<std::endl;
    std::cout<<"Tsecond in: "<<Tsecond_in_combined<<std::endl;
    std::cout<<"phi_g (arg): "<<phi_g<<"\n";
    std::cout<<"phi_g (comp): "<<gas_del.equivalenceRatio()<<"\n";
    for(size_t k=0; k!=nsp; ++k)
        std::cout<<gas.speciesName(k)<<" "<<YSecond_in[k]<<"\n";
    //=====================================================================

    // setup the second stage
    gas.setState_TPY(Tsecond_in_combined, pressure, YSecond_in.data());
    gas.getMoleFractions(x.data());
    rho_in = gas.density();
    gas.getMassFractions(&yin[0]);
    gas.equilibrate("HP");
    Tb = gas.temperature();
    gas.getMassFractions(&yout[0]);
    gas.setState_TP(Tb, pressure);
    gas.getMassFractions(&yout[0]);
    rho_out = gas.density();
    Tad = gas.temperature();
    gas.setState_TPY(Tsecond_in_combined, pressure, YSecond_in.data());



    Cantera::PorousFlow flow_second(foam_stack_secondStage, &gas);
    z.resize(nInitialPoints+2);
    z[0]=0;
    dz = initialFlameProfileThickness / (nInitialPoints - 1);
    expected_flame_location = 0.5*initialFlameProfileThickness+0.005;
    for (std::size_t iz = 0; iz < nInitialPoints; iz++)
        z[iz + 1] = iz*dz + expected_flame_location - initialFlameProfileThickness*0.5;
    z[nInitialPoints+1] = totalLength_second;
    flow_second.setupGrid(z.size(), &z[0]);
    std::unique_ptr<Cantera::Transport> trmix_second(Cantera::newTransportMgr("Mix", sol.thermo().get()));
    flow_second.setTransport(*trmix_second);
    flow_second.setKinetics(*sol.kinetics());
    flow_second.setPressure(pressure);
    flow_second.set_Tamb(Tamb);
    flow_second.set_Tsolid_inlet(Tsfirst_out);
    flow_second.setSolid(!deactivateSolid);
    flow_second.setSolidAxialRadiation(radiation);
    flow_second.setSolidRadialRadiation(radiation);
    flow_second.setAxisymmetricFlow();
    Cantera::InletWithSolid1D inlet_second;
    inlet_second.setFoams(&flow_second.foams);
    inlet_second.setMoleFractions(x.data());
    uin = mdot_g/(rho_in*foam_stack_secondStage.get_local_foam_properties(0.0, Tsecond_in_combined).porosity);
    inlet_second.setMdot(mdot_g);
    inlet_second.setTemperature(Tfirst_out);
    inlet_second.setTsolid(Tsfirst_out);
    Cantera::OutletWithSolid1D outlet_second;
    outlet_second.setFoams(&flow_second.foams);
    outlet_second.setTamb(Tamb);
    std::vector<Cantera::Domain1D*> domains_second { &inlet_second, &flow_second, &outlet_second };
    Cantera::Sim1D flame_second(domains_second);
    flame_second.setGridMin(1,1e-6);
    locs={0.0, (expected_flame_location-0.5*initialFlameProfileThickness)/totalLength_second,
        (expected_flame_location+0.5*initialFlameProfileThickness)/totalLength_second, 1.0};
    porosity_out = foam_stack_secondStage.get_local_foam_properties(totalLength_second, Tad).porosity;
    uout = mdot_g/(rho_out*porosity_out);
    value={uin, uin, uout, uout};
    flame_second.setInitialGuess("velocity",locs,value);
    value = {Tfirst_out, Tfirst_out, Tad, Tad};
    flame_second.setInitialGuess("T",locs,value);
    for (size_t i=0; i<nsp; i++)
    {
        value = {yin[i], yin[i], yout[i], yout[i]};
        flame_second.setInitialGuess(gas.speciesName(i),locs,value);
    }
    value = {Tsecond_in_combined, Tsecond_in_combined, Tad, Tad};
    flame_second.setInitialGuess("Tsolid",locs,value);
    inlet_second.setMoleFractions(x.data());
    inlet_second.setMdot(mdot_g);
    inlet_second.setTemperature(Tsecond_in_combined);
    inlet_second.setTsolid(Tsfirst_out);
    initial_Ts_x = {0, 0.001, 0.002, totalLength_second};
    initial_Ts = {Tsfirst_out, Tsfirst_out, Tad, Tad};
    flow_second.set_initial_guess_Ts(initial_Ts_x, initial_Ts);
    flow_second.set_relaxation(relax);
    flow_second.set_fixedTemperature(true);
    flow_second.solveEnergyEqn();


    auto setTx_second = [&flow_second,&flame_second,flowdomain](Cantera::vector_fp& x, Cantera::vector_fp& T, Cantera::vector_fp& Ts)
    {
        T.resize(flow_second.nPoints());
        Ts.resize(flow_second.nPoints());
        x.resize(flow_second.nPoints());
        for (std::size_t i=0; i!=flow_second.nPoints(); ++i)
        {
            x[i] = flow_second.grid(i);
            Ts[i] = flame_second.value(flowdomain,flow_second.componentIndex("Tsolid"),i);
            T[i] = flame_second.value(flowdomain,flow_second.componentIndex("T"),i);
        }
    };

    // solve the second stage
    {
        flow_second.setTransientTolerances(3e-2, 1e-6);
        flow_second.setSteadyTolerances(3e-2, 1e-6);
        flame_second.setRefineCriteria(flowdomain, 3, 0.5, 0.5, -0.05);
        flame_second.solve(loglevel,true);
        double ratio = 2.0;
        double slope = 0.15;
        double curve = 0.25;
        double prune = 0.01;
        flame_second.setRefineCriteria(flowdomain, ratio, slope, curve, prune);

        flame_second.solve(loglevel,true);
        flow_second.setTransientTolerances(1e-4, 1e-9);
        flow_second.setSteadyTolerances(1e-4, 1e-9);

        flame_second.solve(loglevel,true);

        flow_second.setTransientTolerances(1e-7, 1e-13);
        flow_second.setSteadyTolerances(1e-7, 1e-13);

        Cantera::vector_fp gasT_prev;
        Cantera::vector_fp gasTs_prev;
        Cantera::vector_fp gasx_prev;
        Cantera::vector_fp gasT_curr;
        Cantera::vector_fp gasTs_curr;
        Cantera::vector_fp gasx_curr;
        Cantera::vector_fp T_tmp;

        setTx_second(gasx_prev, gasT_prev, gasTs_prev);
        for (int i=0; i!=iterations; ++i)
        {
            flow_second.set_firstSolve(true);
            flame_second.solve(loglevel,true);

            setTx_second(gasx_curr, gasT_curr, gasTs_curr);
            std::cout<<"Coupling gas and solid: "<<i+1<<"/"<<iterations<<'\n';

            // interpolate old results onto grid of new solution
            interpolate_T(gasx_curr, gasx_prev, gasT_prev, T_tmp);
            double L2gas = L2(gasT_curr, T_tmp);
            std::cout<<"initial solution stage 2: L2 Tg: "<<L2gas<<'\n';
            interpolate_T(gasx_curr, gasx_prev, gasTs_prev, T_tmp);
            double L2solid = L2(gasTs_curr, T_tmp);
            std::cout<<"initial solution stage 2: L2 Ts: "<<L2solid<<'\n';

            if (L2gas < residual_second && L2solid < residual_second)
            {
                std::cout<<'\n'<<"Phases converged!"<<'\n';
                break;
            }
            gasx_curr.swap(gasx_prev);
            gasT_curr.swap(gasT_prev);
            gasTs_curr.swap(gasTs_prev);
        }

        std::cout<<"success"<<std::endl;

        for(size_t i=0;i!=flow_second.nPoints();++i)
        {
            std::cout<<flow_second.grid(i)<<" "<<flame_second.value(flowdomain,flow_second.componentIndex("T"),i)<<" "<<flame_second.value(flowdomain,flow_second.componentIndex("Tsolid"),i)<<"\n";
        }

        /*
        std::cout<<"Inlet conditions for second stage:"<<std::endl;
        std::cout<<"Tfirst out: "<<Tfirst_out<<std::endl;
        std::cout<<"Tsecond in: "<<Tsecond_in_combined<<std::endl;
        std::cout<<"phi_g (arg): "<<phi_g<<"\n";
        std::cout<<"phi_g (comp): "<<gas_del.equivalenceRatio()<<"\n";
        for(size_t k=0; k!=nsp; ++k)
            std::cout<<gas.speciesName(k)<<" "<<YSecond_in[k]<<"\n";
        std::cout<<"Tad "<<Tad<<std::endl;
        */
    }

    // at this point, we have converged solutions for the first stage, and the second stage with correct inlet conditions

    Cantera::vector_fp xc;
    Cantera::vector_fp prevMesh;
    Cantera::vector_fp TsFixed;

    Cantera::vector_fp gasTs_prev;
    Cantera::vector_fp gasT_prev;
    Cantera::vector_fp gasx_prev;
    Cantera::vector_fp gasTs_curr;
    Cantera::vector_fp gasT_curr;
    Cantera::vector_fp gasx_curr;
    Cantera::vector_fp T_tmp;

    Cantera::vector_fp gasTs_prev_second;
    Cantera::vector_fp gasT_prev_second;
    Cantera::vector_fp gasx_prev_second;
    Cantera::vector_fp gasTs_curr_second;
    Cantera::vector_fp gasT_curr_second;
    Cantera::vector_fp gasx_curr_second;

    setTx_second(gasx_prev_second, gasT_prev_second, gasTs_prev_second);
    setTx(gasx_prev, gasT_prev, gasTs_prev);

    for(int it = 0; it != 1000; ++it)
    {

        solveTsolid(flow, flow_second, flame, flame_second, prevMesh, xc, TsFixed); // solve the solid temperature equation for both domains coupled

        //fix the new solid temperatures
        flow.xc.resize(flow.nPoints());
        flow.TsFixed.resize(flow.nPoints());
        for(size_t i = 0; i!=flow.nPoints(); ++i)
        {
            flow.xc[i] = xc[i];
            flow.TsFixed[i] = TsFixed[i];
        }

        flow_second.xc.resize(flow_second.nPoints());
        flow_second.TsFixed.resize(flow_second.nPoints());
        for(size_t i = 0; i!=flow_second.nPoints(); ++i)
        {
            flow_second.xc[i] = xc[i+flow.nPoints()-1];
            flow_second.TsFixed[i] = TsFixed.at(i+flow.nPoints()-1);
        }
        flow_second.xc[0] = 0.;
        // do not solve the solid temperature. use the fixed values
        flow.set_externTsolidEveryIteration(false);
        flow_second.set_externTsolidEveryIteration(false);

        //solve first stage
        flame.solve(loglevel,true);

        Cantera::vector_fp Yfirst_out(nsp);
        for(size_t k=0; k!=nsp; ++k)
            Yfirst_out[k] = flame.value(flowdomain,flow.componentIndex(gas.speciesName(k)),flow.nPoints()-1);
        double Tfirst_out = flame.value(flowdomain,flow.componentIndex("T"),flow.nPoints()-1);
        gas_del.setState_TPY(Tfirst_out,p0,Yfirst_out.data());
        double hfirst_out = gas_del.enthalpy_mass();
        double Tsfirst_out = flame.value(flowdomain,flow.componentIndex("Tsolid"),flow.nPoints()-1);

        Cantera::vector_fp YSecond_in(nsp);
        for(size_t k=0; k!=nsp; ++k)
            YSecond_in[k] = mdot_p * Yfirst_out[k] + mdot_s * Y_second[k];
        gas_del.setState_TPY(Tsecond_in,p0,YSecond_in.data());
        gas_del.getMassFractions(YSecond_in.data());
        gas_del.getMoleFractions(x.data());
        double hsecond_in_combined = mdot_p/mdot_g * hfirst_out + mdot_s/mdot_g*hsecond_in;
        gas_del.setState_HP(hsecond_in_combined, p0);
        double Tsecond_in_combined = gas_del.temperature();

        //update inlet of second stage
        flow_second.set_Tsolid_inlet(Tsfirst_out);
        inlet_second.setMoleFractions(x.data());
        inlet_second.setTemperature(Tsecond_in_combined);
        inlet_second.setTsolid(Tsfirst_out);

        flow_second.TsFixed[0] = Tsfirst_out;

        //solve second stage with updated inlet
        flame_second.solve(loglevel,true);



        // evaluate temperature residuals on each domain
        setTx(gasx_curr, gasT_curr, gasTs_curr);
        interpolate_T(gasx_curr, gasx_prev, gasT_prev, T_tmp);
        double L2gas= L2(gasT_curr, T_tmp);
        interpolate_T(gasx_curr, gasx_prev, gasTs_prev, T_tmp);
        double L2solid= L2(gasTs_curr, T_tmp);

        setTx_second(gasx_curr_second, gasT_curr_second, gasTs_curr_second);
        interpolate_T(gasx_curr_second, gasx_prev_second, gasT_prev_second, T_tmp);
        double L2gas_second = L2(gasT_curr_second, T_tmp);
        interpolate_T(gasx_curr_second, gasx_prev_second, gasTs_prev_second, T_tmp);
        double L2solid_second = L2(gasTs_curr_second, T_tmp);
        std::cout<<"L2 1st stage gas:   "<<L2gas<<'\n';
        std::cout<<"L2 1st stage solid: "<<L2solid<<'\n';
        std::cout<<"L2 2nd stage gas:   "<<L2gas_second<<'\n';
        std::cout<<"L2 2nd stage solid: "<<L2solid_second<<'\n';
        gasx_curr.swap(gasx_prev);
        gasT_curr.swap(gasT_prev);
        gasTs_curr.swap(gasTs_prev);
        gasx_curr_second.swap(gasx_prev_second);
        gasT_curr_second.swap(gasT_prev_second);
        gasTs_curr_second.swap(gasTs_prev_second);

        if (L2gas<final_residual && L2solid<final_residual && L2gas_second<final_residual && L2solid_second<final_residual)
            break;
    }

    std::ofstream write("table.txt");
    write<<"#x u rho mdot phi T Ts ";
    for(size_t k=0; k!=gas.nSpecies(); ++k)
        write << gas.speciesName(k)<<" ";
    write<<'\n';

    vector_fp Ytmp(nsp);
    for(size_t i=0;i!=flow.nPoints()-1; ++i)
    {
        for(size_t k=0;k!=nsp;++k)
            Ytmp[k] = flame.value(flowdomain,flow.componentIndex(gas.speciesName(k)),i);
        double Ttmp = flame.value(flowdomain,flow.componentIndex("T"),i);
        gas.setState_TPY(Ttmp,p0,Ytmp.data());

        write << flow.grid(i) << " ";
        write << flame.value(flowdomain,flow.componentIndex("velocity"),i) << " ";
        write << gas.density() << " ";
        write << flame.value(flowdomain,flow.componentIndex("velocity"),i) * gas.density() * foam_stack.get_local_foam_properties(flow.grid(i), Ttmp).porosity<<" ";
        write << gas.equivalenceRatio() << " ";
        write << Ttmp<<" ";
        write << flame.value(flowdomain,flow.componentIndex("Tsolid"),i)<<" ";
        for(size_t k=0;k!=nsp;++k)
            write << Ytmp[k] << " ";
        write<<'\n';
    }

    for(size_t i=0;i!=flow_second.nPoints(); ++i)
    {
        for(size_t k=0;k!=nsp;++k)
            Ytmp[k] = flame_second.value(flowdomain,flow_second.componentIndex(gas.speciesName(k)),i);
        double Ttmp = flame_second.value(flowdomain,flow_second.componentIndex("T"),i);
        gas.setState_TPY(Ttmp,p0,Ytmp.data());

        write << flow_second.grid(i)+totalLength << " ";
        write << flame_second.value(flowdomain,flow_second.componentIndex("velocity"),i) << " ";
        write << gas.density() << " ";
        write << flame_second.value(flowdomain,flow_second.componentIndex("velocity"),i) * gas.density() * foam_stack_secondStage.get_local_foam_properties(flow_second.grid(i), Ttmp).porosity<<" ";
        write << gas.equivalenceRatio() << " ";
        write << Ttmp<<" ";
        write << flame_second.value(flowdomain,flow_second.componentIndex("Tsolid"),i)<<" ";
        for(size_t k=0;k!=nsp;++k)
            write << Ytmp[k] << " ";
        write<<'\n';
    }
    write.close();
    std::exit(EXIT_SUCCESS);
}

