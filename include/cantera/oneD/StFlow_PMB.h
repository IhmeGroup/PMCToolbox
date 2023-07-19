//! @file StFlow_PMB.h

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

#ifndef CT_STFLOW_PMB_H
#define CT_STFLOW_PMB_H

#include "StFlow.h"
#include "cantera/transport.h"
#include "cantera/numerics/eigen_sparse.h"
#include "cantera/base/Solution.h"
#include "cantera/numerics/Func1.h"
#include <fstream>

namespace Cantera
{

class FoamProperties
{
    public:

        FoamProperties() : gap(false) {}
        FoamProperties(const FoamProperties&) = default;

        bool gap; // determines if this represents a gap between solid matrices or a solid matrix. First and last block cannot be a gap
                  // gap regons interpolate properties between to solid matrices linearly to avoid a jump in properties

        /* matrix properties */
        std::shared_ptr<Cantera::Func1> porosity; // porosity as a function of local x (m)
        std::shared_ptr<Cantera::Func1> heat_conductivity; // solid heat conductivity W/m/K as a function of T (K)
        std::shared_ptr<Cantera::Func1> specific_area; // specific area (1/m) as a function of x (m)
        std::shared_ptr<Cantera::Func1> heatCond_eff_factor; // effectivity factor for solid heat conduction 
        std::shared_ptr<Cantera::Func1> tortuosity_factor; // tortuosity factor for species and heat diffusion
        std::shared_ptr<Cantera::Func1> extinction_coefficient; // for solid radiation modeling
        std::shared_ptr<Cantera::Func1> hydraulic_diameter; // for Reynolds number

        /* material properties */
        std::shared_ptr<Cantera::Solution> solidPhase; // used to compute density and heat capacity of the solid

        double height; // height of the solid matrix section

        /* Radiation properties */
        double emissivity; // only used if doSolidRadialRadiation = true and/or doSolidAxialRadiation = true
        double insulationTransmissivity; // only used if doSolidRadialRadiation = true
        double chemistryFactor; // optionally, multiply chemical reaction rates and heat release rate in this section by this factor (default: 1)

        static FoamProperties createGap(double height)
        {
            FoamProperties prop;
            prop.gap = true;
            prop.height = height;
            return prop;
        }
};


class FoamStack
{

    int get_foam_index(double x)
    {
        for (int i=0; i!=static_cast<int>(foams.size()); ++i)
        {
            const FoamProperties& foam = foams[i];
            if (x <= x_start[i] + foam.height)
                return i;
        }
        return foams.size()-1;
    }

    public:
        FoamStack() = delete;
        FoamStack(const FoamStack&) = delete;
        FoamStack& operator=(const FoamStack&) = delete;
        FoamStack(const std::vector<FoamProperties>& foamList, double diam) :
            foams(foamList), burner_diameter(diam)
        {
            if (foams.empty())
            {
                std::cout<<"ERROR: empty list of foams!"<<std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (foams[0].gap || foams[foams.size()-1].gap)
            {
                std::cout<<"ERROR: first and/or last foam cannot be a gap!"<<std::endl;
                std::exit(EXIT_FAILURE);
            }
            double sum = 0.0;
            for (const auto& e : foams)
            {
                x_start.push_back(sum);
                sum += e.height;
            }
            totalHeight = sum;
        }

        struct LocalFoamProperties
        {
            double density;
            double effective_heat_conductivity;
            double heat_capacity;

            double porosity;

            double tortuosity_factor;
            double scatteringAlbedo;
            double extinction_coefficient;

            double hydraulic_diameter;
            double specific_area;

            double emissivity;
            double insulationTransmissivity;
            double chemistryFactor;
        };


        LocalFoamProperties get_local_foam_properties(double localx, int index, double Ts, double Tmid = -1, int id=-123456)
        {
            if (Ts <= 250.)
            {
                Ts = 300;
            }
            if (Tmid < 0.)
                Tmid = Ts;

            const FoamProperties& foam = foams[index];
            LocalFoamProperties prop;

            auto& thermoPhase = *foam.solidPhase->thermo();
            thermoPhase.setState_TP(Ts, thermoPhase.pressure());

            prop.density               = thermoPhase.density();
            prop.heat_capacity         = thermoPhase.cp_mass();

            prop.porosity              = (*(foam.porosity))(localx);
            prop.tortuosity_factor     = (*(foam.tortuosity_factor))(localx);
            prop.extinction_coefficient = (*(foam.extinction_coefficient))(localx);
            prop.specific_area         = (*(foam.specific_area))(localx);

            prop.scatteringAlbedo      = 0.5*(2.0-foam.emissivity);
            prop.extinction_coefficient = (*(foam.extinction_coefficient))(localx);
            prop.hydraulic_diameter     = (*(foam.hydraulic_diameter))(localx);

            prop.effective_heat_conductivity     = ((*(foam.heatCond_eff_factor))(localx))*((*(foam.heat_conductivity))(Tmid));

            prop.emissivity = foam.emissivity;
            prop.insulationTransmissivity = foam.insulationTransmissivity;
            prop.chemistryFactor = foam.chemistryFactor;
            return prop;

        }

        LocalFoamProperties get_local_foam_properties(double globalx, double Ts, double Tmid = -1.0, int id=-123456)
        {
            if (Tmid < 0.0)
                Tmid = Ts;
            int index = get_foam_index(globalx);
            double localx = globalx - x_start[index];
            const FoamProperties& foam = foams[index];

            LocalFoamProperties prop;
            if (foam.gap)
            {
                auto linInterpol = [](double x, double x0, double x1, double y0, double y1) { return y0 + (x - x0)*(y1-y0)/(x1-x0); };

                const auto& left  = get_local_foam_properties(foams[index-1].height, index-1, Ts, Tmid, id);
                const auto& right = get_local_foam_properties(0.0, index+1, Ts, Tmid, id);

                prop.density                        = linInterpol(localx, 0.0, foam.height, left.density,                       right.density);
                prop.effective_heat_conductivity    = linInterpol(localx, 0.0, foam.height, left.effective_heat_conductivity,   right.effective_heat_conductivity);
                prop.heat_capacity                  = linInterpol(localx, 0.0, foam.height, left.heat_capacity,                 right.heat_capacity);
                prop.porosity                       = linInterpol(localx, 0.0, foam.height, left.porosity,                      right.porosity);
                prop.scatteringAlbedo               = linInterpol(localx, 0.0, foam.height, left.scatteringAlbedo,              right.scatteringAlbedo);
                prop.tortuosity_factor              = linInterpol(localx, 0.0, foam.height, left.tortuosity_factor,             right.tortuosity_factor);
                prop.extinction_coefficient         = linInterpol(localx, 0.0, foam.height, left.extinction_coefficient,        right.extinction_coefficient);
                prop.hydraulic_diameter             = linInterpol(localx, 0.0, foam.height, left.hydraulic_diameter,            right.hydraulic_diameter);
                prop.specific_area                  = linInterpol(localx, 0.0, foam.height, left.specific_area,                 right.specific_area);
                prop.emissivity                     = linInterpol(localx, 0.0, foam.height, left.emissivity,                    right.emissivity);
                prop.insulationTransmissivity       = linInterpol(localx, 0.0, foam.height, left.insulationTransmissivity,      right.insulationTransmissivity);
                prop.chemistryFactor                = linInterpol(localx, 0.0, foam.height, left.chemistryFactor,               right.chemistryFactor);
            }
            else
            {
                prop = get_local_foam_properties(localx, index, Ts, Tmid, id);
            }
            return prop;
        }

    public:
        std::vector<FoamProperties> foams; // foam stack from bottom to top
        std::vector<double> x_start;
        double totalHeight;
        double burner_diameter;
};


class PorousFlow : public StFlow
{
public:

    void setSolid(bool b){ m_doSolid=b; }
    bool doSolid(){ return m_doSolid; }

    void setSolidAxialRadiation(bool b){ m_doSolidAxialRadiation=b; }
    bool doSolidAxialRadiation(){ return m_doSolidAxialRadiation; }

    void setSolidRadialRadiation(bool b){ m_doSolidRadialRadiation=b; }
    bool doSolidRadialRadiation(){ return m_doSolidRadialRadiation; }

    void set_Tamb(double T){ m_T_amb = T; }
    double get_Tamb(){ return m_T_amb; }
    void set_Tsolid_inlet(double T){ m_Tsolid_inlet = T; }
    double set_Tsolid_inlet() const { return m_Tsolid_inlet; }


    void setSolidPseudoTimeStepping(bool b){ m_solidPseudoTimeStepping=b; }
    void setTsolidZeroGradient(bool b){ m_TsolidZeroGradient=b; }

    doublereal divHeatFluxPorosity(const doublereal* x, size_t j) const {
        doublereal c1 = (m_porosity[j-1]/m_tortuosity_factor[j-1])*m_tcon[j-1]*(T(x,j) - T(x,j-1));
        doublereal c2 = (m_porosity[j]/m_tortuosity_factor[j])  *m_tcon[j]  *(T(x,j+1) - T(x,j));
        return -2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }
    doublereal divHeatFluxSolid(const doublereal* x, size_t j) const {
        doublereal c1 = m_effectiveSolidHeatConductivity[j-1]*(Tsolid(x,j)   - Tsolid(x,j-1));
        doublereal c2 = m_effectiveSolidHeatConductivity[j]  *(Tsolid(x,j+1) - Tsolid(x,j));
        return -2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    PorousFlow(FoamStack& foamList, ThermoPhase* ph = 0, size_t nsp = 1, size_t points = 1);

    void resize(size_t components, size_t points)
    {
        StFlow::resize(components, points);
        m_do_solid_temperature.resize(m_points,false);

        m_divq.resize(m_points, 0.0);
        m_divDiffusiveSolidHeatFlux.resize(m_points, 0.0);
        m_rad_radial.resize(m_points, 0.0);
        m_porosity.resize(m_points, 0.0);
        m_solidDensity.resize(m_points, 0.0);
        m_effectiveSolidHeatConductivity.resize(m_points, 0.0);
        m_solidHeatCapacity.resize(m_points, 0.0);
        m_tortuosity_factor.resize(m_points, 0.0);
        m_hydraulicDiameter.resize(m_points, 0.0);
        m_specificArea.resize(m_points, 0.0);
        m_extinctionCoefficient.resize(m_points, 0.0);
        m_scatteringAlbedo.resize(m_points, 0.0);
        m_emissivity.resize(m_points, 0.0);
        m_hv.resize(m_points, 0.0);
        m_insulationTransmissivity.resize(m_points, 0.0);
        m_chemistryFactor.resize(m_points, 0.0);
    }

    bool componentActive(size_t n) const
    {
        if (n == c_offset_Tsolid)
            return  true;
        else
            return StFlow::componentActive(n);
    }
    void updateTransport(double* c, size_t j0, size_t j1)
    {
        StFlow::updateTransport(c,j0,j1);
    }

    void updateDiffFluxes(const doublereal* x, size_t j0, size_t j1)
    {
        if (m_do_multicomponent) {
            for (size_t j = j0; j < j1; j++) {
                double dz = z(j+1) - z(j);
                for (size_t k = 0; k < m_nsp; k++) {
                    doublereal sum = 0.0;
                    for (size_t m = 0; m < m_nsp; m++) {
                        sum += m_wt[m] * m_multidiff[mindex(k,m,j)] * (X(x,m,j+1)-X(x,m,j));
                    }
                    m_flux(k,j) = sum * m_diff[k+j*m_nsp] / dz;
                }
            }
        } else {
            for (size_t j = j0; j < j1; j++) {
                double sum = 0.0;
                double wtm = m_wtm[j];
                double rho = density(j);
                double dz = z(j+1) - z(j);
                for (size_t k = 0; k < m_nsp; k++) {
                    m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm);
                    m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
                    sum -= m_flux(k,j);
                }
                // correction flux to insure that \sum_k Y_k V_k = 0.
                for (size_t k = 0; k < m_nsp; k++) {
                    m_flux(k,j) += sum*Y(x,k,j);
                }
            }
        }

        if (m_do_soret) {
            for (size_t m = j0; m < j1; m++) {
                double gradlogT = 2.0 * (T(x,m+1) - T(x,m)) /
                                  ((T(x,m+1) + T(x,m)) * (z(m+1) - z(m)));
                for (size_t k = 0; k < m_nsp; k++) {
                    m_flux(k,m) -= m_dthermal(k,m)*gradlogT;
                }
            }
        }
    }



    void evalRightBoundary(double* x, double* rsd, int* diag, double rdt)
    {
        size_t j = m_points - 1;

        // the boundary object connected to the right of this one may modify or
        // replace these equations. The default boundary conditions are zero u, V,
        // and T, and zero diffusive flux for all species.

        rsd[index(c_offset_V,j)] = V(x,j);
        doublereal sum = 0.0;
        rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
        diag[index(c_offset_L, j)] = 0;
        for (size_t k = 0; k < m_nsp; k++) {
            sum += Y(x,k,j);
            rsd[index(k+c_offset_Y,j)] = m_porosity[j]*((1.0/m_tortuosity_factor[j])*m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j));
        }
        rsd[index(c_offset_Y + rightExcessSpecies(), j)] = 1.0 - sum;
        diag[index(c_offset_Y + rightExcessSpecies(), j)] = 0;
        if (domainType() == cAxisymmetricStagnationFlow) {
            rsd[index(c_offset_U,j)] = m_porosity[j]*rho_u(x,j);
            if (m_do_energy[j]) {
                rsd[index(c_offset_T,j)] = T(x,j);
            } else {
                rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
            }
        } else if (domainType() == cFreeFlow) {
            rsd[index(c_offset_U,j)] = m_porosity[j]*rho_u(x,j) - m_porosity[j-1]*rho_u(x,j-1);
            rsd[index(c_offset_T,j)] = T(x,j) - T(x,j-1);
        }

        if (m_TsolidZeroGradient)
        {
            rsd[index(c_offset_Tsolid,j)] = Tsolid(x,j);
        }
        else
        {
            rsd[index(c_offset_Tsolid,j)] = (Tsolid(x,j) - Tsolid(x,j-1))*0.5*(m_effectiveSolidHeatConductivity[j]+m_effectiveSolidHeatConductivity[j-1])/(z(j)-z(j-1)) + m_emissivity[j]*StefanBoltz*(1.-m_porosity[j])*(std::pow(Tsolid(x,j),4.0) - std::pow(m_T_amb,4.0));
        }
    }

    void evalResidual(double* x, double* rsd, int* diag,
                      double rdt, size_t jmin, size_t jmax)
    {
        bool recalcMatrixGlobal = false;
        if (m_doSolidAxialRadiation && m_doSolid)
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
                    if (std::abs(z(i) - prevMesh[i]) > 1e-8)
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
                    prevMesh[i] = z(i);
            }

            recalcMatrixGlobal = recalcMatrix;

            constexpr int non_zeros_per_col = 3;

            if (recalcMatrix)
            {
                mat = Eigen::SparseMatrix<double,Eigen::RowMajor>(2*m_points,2*m_points);
                mat.reserve(Eigen::Matrix<int, Eigen::Dynamic, 1>::Constant(2*m_points, non_zeros_per_col));

                mat.insert(0,0) = 1.0;
                for (int i=1; i!=static_cast<int>(m_points); ++i)
                {
                    mat.insert(i, i-1) = -1.0;
                    mat.insert(i, i) = 1.0+(z(i)-z(i-1))*m_extinctionCoefficient[i]*(2.0-m_scatteringAlbedo[i]);
                    mat.insert(i, i+m_points) = -(z(i)-z(i-1))*m_extinctionCoefficient[i]*m_scatteringAlbedo[i];
                }
                for (int j=m_points; j!=2*static_cast<int>(m_points)-1; ++j)
                {
                    int i = j-m_points;
                    mat.insert(j, j-m_points) = -(z(i+1)-z(i))*m_extinctionCoefficient[i]*m_scatteringAlbedo[i];
                    mat.insert(j, j) = 1.0+(z(i+1)-z(i))*m_extinctionCoefficient[i]*(2.0-m_scatteringAlbedo[i]);
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
                radiation_rhs(i) = (z(i)-z(i-1))*2.0*m_extinctionCoefficient[i]*(1.0-m_scatteringAlbedo[i])*StefanBoltz*std::pow(Tsolid(x,i), 4.0);
            for (int j=m_points; j!=2*static_cast<int>(m_points)-1; ++j)
            {
                int i = j-m_points;
                radiation_rhs(j) = (z(i+1)-z(i))*2.0*m_extinctionCoefficient[i]*(1.0-m_scatteringAlbedo[i])*StefanBoltz*std::pow(Tsolid(x,i), 4.0);
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
                m_divq[i] = 2*m_extinctionCoefficient[i]*(1.0-m_scatteringAlbedo[i])*(2.0*StefanBoltz*std::pow(Tsolid(x,i), 4.0) - (q[i]+q[i+m_points]));

        }


        if (m_doSolidRadialRadiation)
        {
            for (int i=0; i!=static_cast<int>(m_points); ++i)
                m_rad_radial[i] = 4.0*StefanBoltz*m_emissivity[i]*m_insulationTransmissivity[i]*(std::pow(Tsolid(x,i), 4.0) - std::pow(m_T_amb, 4.0))/(foams.burner_diameter);
        }


        if (m_externTsolidEveryIteration)
        {
            if (m_points > 35 && (recalcMatrixGlobal || m_firstSolve))
            {
                m_firstSolve = false;
                Tgprev.resize(m_points);
                Tsprev.resize(m_points);
                for(size_t i=0;i!=m_points; ++i)
                {
                    Tgprev[i] = T(x,i);
                    Tsprev[i] = Tsolid(x,i);
                }
                solveTsolid(Tgprev,Tsprev);
            }
        }

        for (size_t j = jmin; j <= jmax; j++) {
            m_divDiffusiveSolidHeatFlux[j] = divHeatFluxSolid(x,j);
            //----------------------------------------------
            //         left boundary
            //----------------------------------------------
            if (j == 0) {
                // Continuity. This propagates information right-to-left, since
                // rho_u at point 0 is dependent on rho_u at point 1, but not on
                // mdot from the inlet.
                rsd[index(c_offset_U,0)] =
                    -(m_porosity[1]*rho_u(x,1) - m_porosity[0]*rho_u(x,0))/m_dz[0]
                    -(m_porosity[1]*density(1)*V(x,1) - m_porosity[0]*density(0)*V(x,0));

                // the inlet (or other) object connected to this one will modify
                // these equations by subtracting its values for V, T, and mdot. As
                // a result, these residual equations will force the solution
                // variables to the values for the boundary object

                if (m_doSolid)
                {
                    rsd[index(c_offset_Tsolid, j)] = Tsolid(x,0);
                    diag[index(c_offset_Tsolid, j)] = 1;
                }
                else
                {
                    rsd[index(c_offset_Tsolid, j)] = Tsolid(x,0);
                    diag[index(c_offset_Tsolid, j)] = 0;
                }

                rsd[index(c_offset_V,0)] = V(x,0);
                if (doEnergy(0)) {
                    rsd[index(c_offset_T,0)] = T(x,0);
                } else {
                    rsd[index(c_offset_T,0)] = T(x,0) - T_fixed(0);
                }
                rsd[index(c_offset_L,0)] = -m_porosity[0]*rho_u(x,0);

                // The default boundary condition for species is zero flux. However,
                // the boundary object may modify this.
                double sum = 0.0;
                for (size_t k = 0; k < m_nsp; k++) {
                    sum += Y(x,k,0);
                    rsd[index(c_offset_Y + k, 0)] =
                        -((m_porosity[0]/m_tortuosity_factor[0])*m_flux(k,0) + m_porosity[0]*rho_u(x,0)* Y(x,k,0));
                }
                rsd[index(c_offset_Y + leftExcessSpecies(), 0)] = 1.0 - sum;

                // set residual of poisson's equ to zero
                rsd[index(c_offset_E, 0)] = x[index(c_offset_E, j)];
            } else if (j == m_points - 1) {
                evalRightBoundary(x, rsd, diag, rdt);
                // set residual of poisson's equ to zero
                rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];

                if (m_doSolid)
                {
                    if (m_TsolidZeroGradient)
                    {
                        rsd[index(c_offset_Tsolid,j)] = Tsolid(x,j);
                    }
                    else
                    {
                        rsd[index(c_offset_Tsolid,j)] = (Tsolid(x,j) - Tsolid(x,j-1))*0.5*(m_effectiveSolidHeatConductivity[j]+m_effectiveSolidHeatConductivity[j-1])/(z(j)-z(j-1)) + m_emissivity[j]*StefanBoltz*(1.-m_porosity[j])*(std::pow(Tsolid(x,j),4.0) - std::pow(m_T_amb,4.0));
                    }
                    diag[index(c_offset_Tsolid, j)] = 1;
                }
                else
                {
                    rsd[index(c_offset_Tsolid,j)] = Tsolid(x,j);
                    diag[index(c_offset_Tsolid, j)] = 0; //?
                }

            } else { // interior points

                evalContinuity(j, x, rsd, diag, rdt);
                // set residual of poisson's equ to zero
                rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];

                //------------------------------------------------
                //    Radial momentum equation
                //
                //    \rho dV/dt + \rho u dV/dz + \rho V^2
                //       = d(\mu dV/dz)/dz - lambda
                //-------------------------------------------------
                rsd[index(c_offset_V,j)]
                    = (shear(x,j) - lambda(x,j) - rho_u(x,j)*dVdz(x,j)
                            - m_rho[j]*V(x,j)*V(x,j))/m_rho[j]
                    - rdt*(V(x,j) - V_prev(j));
                diag[index(c_offset_V, j)] = 1;

                //-------------------------------------------------
                //    Species equations
                //
                //   \rho dY_k/dt + \rho u dY_k/dz + dJ_k/dz
                //   = M_k\omega_k
                //-------------------------------------------------
                getWdot(x,j);
                for (size_t k = 0; k < m_nsp; k++) {
                    double convec = m_porosity[j]*rho_u(x,j)*dYdz(x,k,j);
                    double diffus = 2.0*((m_porosity[j]/m_tortuosity_factor[j])*m_flux(k,j) - (m_porosity[j-1]/m_tortuosity_factor[j-1])*m_flux(k,j-1))
                        / (z(j+1) - z(j-1));
                    rsd[index(c_offset_Y + k, j)]
                        = (m_porosity[j]*m_wt[k]*(m_chemistryFactor[j]*wdot(k,j))
                                - convec - diffus)/(m_porosity[j]*m_rho[j])
                        - rdt*(Y(x,k,j) - Y_prev(k,j));
                    diag[index(c_offset_Y + k, j)] = 1;
                }

                //-----------------------------------------------
                //    energy equation
                //
                //    \rho c_p dT/dt + \rho c_p u dT/dz
                //    = d(k dT/dz)/dz
                //      - sum_k(\omega_k h_k_ref)
                //      - sum_k(J_k c_p_k / M_k) dT/dz
                //-----------------------------------------------
                if (m_do_energy[j]) {
                    setGas(x,j);

                    // heat release term
                    const vector_fp& h_RT = m_thermo->enthalpy_RT_ref();
                    const vector_fp& cp_R = m_thermo->cp_R_ref();
                    double sum = 0.0;
                    double sum2 = 0.0;
                    for (size_t k = 0; k < m_nsp; k++) {
                        double flxk = 0.5*((m_porosity[j-1]/m_tortuosity_factor[j-1])*m_flux(k,j-1) + (m_porosity[j]/m_tortuosity_factor[j])*m_flux(k,j));
                        sum += m_porosity[j]*m_chemistryFactor[j]*wdot(k,j)*h_RT[k];
                        sum2 += flxk*cp_R[k]/m_wt[k];
                    }
                    sum *= GasConstant * T(x,j);
                    double dtdzj = dTdz(x,j);
                    sum2 *= GasConstant * dtdzj;

                    rsd[index(c_offset_T, j)] = - m_porosity[j]*m_cp[j]*rho_u(x,j)*dtdzj
                        - divHeatFluxPorosity(x,j) - sum - sum2;
                    rsd[index(c_offset_T, j)] /= (m_porosity[j]*m_rho[j]*m_cp[j]);
                    rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j));
                    rsd[index(c_offset_T, j)] -= (m_qdotRadiation[j] / (m_rho[j] * m_cp[j]));
                    diag[index(c_offset_T, j)] = 1;

                } else {
                    // residual equations if the energy equation is disabled
                    rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
                    diag[index(c_offset_T, j)] = 0;
                }

                rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
                diag[index(c_offset_L, j)] = 0;

                double gasEnergyExchange = -m_hv[j] * (T(x,j) - Tsolid(x,j));
                if (m_doSolid)
                {
                    rsd[index(c_offset_T, j)] += gasEnergyExchange/(m_porosity[j]*m_rho[j]*m_cp[j]);
                }

                if (!fixedTemperature)
                {
                    if (m_doSolid)
                    {
                        rsd[index(c_offset_Tsolid, j)] = -divHeatFluxSolid(x,j) - gasEnergyExchange - m_divq[j] - m_rad_radial[j];
                        if (m_solidPseudoTimeStepping)
                        {
                            rsd[index(c_offset_Tsolid, j)] /= (1.0-m_porosity[j])*m_solidDensity[j]*m_solidHeatCapacity[j];
                            rsd[index(c_offset_Tsolid, j)] -= rdt*(Tsolid(x,j) - Tsolid_prev(j));
                        }
                        diag[index(c_offset_Tsolid, j)] = 1;
                    }
                    else
                    {
                        rsd[index(c_offset_Tsolid, j)] = 0;
                        diag[index(c_offset_Tsolid, j)] = 0;

                    }
                }
            }
        }

        if (fixedTemperature)
        {
            for (size_t j = jmin; j <= jmax; j++)
            {
                double Tinterpolated = 300.;
                if (z(j) <= xc[0])
                    Tinterpolated = TsFixed[0];
                else if (z(j) >= xc[xc.size()-1])
                    Tinterpolated = TsFixed[TsFixed.size()-1];
                else
                {
                    auto ind = std::lower_bound(xc.begin(), xc.end(), z(j)) - xc.begin();
                    Tinterpolated = TsFixed[ind-1] + (z(j)-xc[ind-1])*(TsFixed[ind]-TsFixed[ind-1])/(xc[ind]-xc[ind-1]);
                }
                rsd[index(c_offset_Tsolid, j)] = Tsolid(x,j) - Tinterpolated;
                diag[index(c_offset_Tsolid, j)] = 0; //?
            }
        }
    }

    void solve_solid_temperature(size_t j=npos);

    void finalize(const double* x)
    {
        StFlow::_finalize(x);

        bool p = m_do_solid_temperature[0];
        if (p) {
            solve_solid_temperature();
        }
    }

    double Tsolid(const double* x, size_t j) const
    {
        return x[index(c_offset_Tsolid, j)];
    }

    double Tsolid_prev(size_t j) const {
        return prevSoln(c_offset_Tsolid, j);
    }

    bool doSolidTemperature(size_t j) {
        return m_do_solid_temperature[j];
    }

    void set_externTsolidEveryIteration(bool b){ m_externTsolidEveryIteration=true; }
    void set_relaxation(double r){m_relax=r;}

    virtual void updateSolidPhaseProperties(double* x, size_t jmin, size_t jmax)
    {
        // properties are computed for grid points from j0 to j1
        size_t j0 = std::max<size_t>(jmin, 1) - 1;
        size_t j1 = std::min(jmax+1,m_points-1);

        // update properties for the solid phase
        for (size_t j = j0; j <= j1; j++)
        {

            double Tmid = (j == m_points-1) ? Tsolid(x,j) : 0.5*(Tsolid(x,j)+Tsolid(x,j+1));
            const auto& solidProperties = foams.get_local_foam_properties(z(j),Tsolid(x,j), Tmid, j);

            m_solidDensity[j]          = solidProperties.density;
            m_solidHeatCapacity[j]     = solidProperties.heat_capacity;
            m_porosity[j]              = solidProperties.porosity;
            m_tortuosity_factor[j]     = solidProperties.tortuosity_factor;
            m_hydraulicDiameter[j]     = solidProperties.hydraulic_diameter;
            m_specificArea[j]          = solidProperties.specific_area;
            m_chemistryFactor[j]       = solidProperties.chemistryFactor;
            m_insulationTransmissivity[j] = solidProperties.insulationTransmissivity;
            m_emissivity[j]               = solidProperties.emissivity;

            setGas(x,j);
            double visc = m_trans->viscosity();
            double cond = m_trans->thermalConductivity();
            double hydraulic_diameter = solidProperties.hydraulic_diameter;
            double Re = std::abs(rho_u(x,j))*hydraulic_diameter/visc;
            double Pr = m_cp[j] * visc / cond;
            double Nu = 3.7*std::pow(Re,0.38)*std::pow(Pr,0.25);
            m_hv[j]                             = Nu*solidProperties.specific_area*cond/hydraulic_diameter;

            m_effectiveSolidHeatConductivity[j] = solidProperties.effective_heat_conductivity;
            m_extinctionCoefficient[j]          = solidProperties.extinction_coefficient;
            m_scatteringAlbedo[j]               = solidProperties.scatteringAlbedo;
        }
    }

    virtual void updateProperties(size_t jg, double* x, size_t jmin, size_t jmax)
    {
        StFlow::updateProperties(jg,x,jmin,jmax);
        updateSolidPhaseProperties(x,jmin,jmax);
    }

void evalContinuity(size_t j, double* x, double* rsd, int* diag, double rdt)
{
    //algebraic constraint
    diag[index(c_offset_U, j)] = 0;
    //----------------------------------------------
    //    Continuity equation
    //
    //    d(\rho u)/dz + 2\rho V = 0
    //----------------------------------------------
    if (domainType() == cAxisymmetricStagnationFlow) {
        // Note that this propagates the mass flow rate information to the left
        // (j+1 -> j) from the value specified at the right boundary. The
        // lambda information propagates in the opposite direction.
        rsd[index(c_offset_U,j)] =
            -(m_porosity[j+1]*rho_u(x,j+1) - m_porosity[j]*rho_u(x,j))/m_dz[j]
            -(density(j+1)*V(x,j+1) + density(j)*V(x,j));
    } else if (domainType() == cFreeFlow) {
        if (grid(j) > m_zfixed) {
            rsd[index(c_offset_U,j)] =
                - (m_porosity[j+1]*rho_u(x,j) - m_porosity[j+1]*rho_u(x,j-1))/m_dz[j-1]
                - (m_porosity[j+1]*density(j-1)*V(x,j-1) + m_porosity[j+1]*density(j)*V(x,j));
        } else if (grid(j) == m_zfixed) {
            if (m_do_energy[j]) {
                rsd[index(c_offset_U,j)] = (T(x,j) - m_tfixed);
            } else {
                rsd[index(c_offset_U,j)] = (m_porosity[j]*rho_u(x,j)
                                            - m_porosity[j]*m_rho[0]*0.3);
            }
        } else if (grid(j) < m_zfixed) {
            rsd[index(c_offset_U,j)] =
                - (m_porosity[j+1]*rho_u(x,j+1) - m_porosity[j]*rho_u(x,j))/m_dz[j]
                - (m_porosity[j+1]*density(j+1)*V(x,j+1) + m_porosity[j+1]*density(j)*V(x,j));
        }
    }
}

void solveEnergyEqn(size_t j=npos);
void fixTemperature(size_t j=npos);

void debugOutput(const std::string& filename, const std::vector<Cantera::vector_fp>& Ys);


public:
void solveTsolid(const std::vector<double>& Tg, const std::vector<double>&Tsolid_prev)
{
    Tsolid_mat = Eigen::SparseMatrix<double,Eigen::RowMajor>(m_points,m_points);
    constexpr int non_zeros_per_col = 3;
    Tsolid_mat.reserve(Eigen::Matrix<int, Eigen::Dynamic, 1>::Constant(m_points, non_zeros_per_col));

    Tsolid_mat.insert(0,0) = 1.0;
    for (int i=1; i!=static_cast<int>(m_points-1); ++i)
    {

        double lambda_mid_left = 0.5*(m_effectiveSolidHeatConductivity[i-1]+m_effectiveSolidHeatConductivity[i]);
        double lambda_mid_right = 0.5*(m_effectiveSolidHeatConductivity[i+1]+m_effectiveSolidHeatConductivity[i]);

        double c1 = 2.0*lambda_mid_left /((z(i+1)-z(i-1))*(z(i)-z(i-1)));
        double c2 = 2.0*lambda_mid_right/((z(i+1)-z(i-1))*(z(i+1)-z(i)));

        Tsolid_mat.insert(i, i-1) = -c1;
        Tsolid_mat.insert(i, i) = m_hv[i] + c1 + c2;
        Tsolid_mat.insert(i, i+1) = -c2;
    }

    if (m_TsolidZeroGradient)
    {
        Tsolid_mat.insert(m_points-1,m_points-1) = 1;
        Tsolid_mat.insert(m_points-1,m_points-2) = -1;
    }
    else
    {
        double lambda_end_mid = 0.5*(m_effectiveSolidHeatConductivity[m_points-1]+m_effectiveSolidHeatConductivity[m_points-2])/(z(m_points-1)-z(m_points-2));
        Tsolid_mat.insert(m_points-1,m_points-1) = lambda_end_mid + m_emissivity[m_points-1]*StefanBoltz*(1.-m_porosity[m_points-1])*Tsolid_prev[m_points-1]*Tsolid_prev[m_points-1]*Tsolid_prev[m_points-1];
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

    if (m_TsolidZeroGradient)
    {
        Tsolid_rhs(m_points-1) = 0;//Tsolid_rhs(m_points-2);
    }
    else
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

    for (int i=0; i!=static_cast<int>(m_points); ++i)
    {
        TsFixed[i] = m_relax*Tsolid_q(i) + (1.-m_relax)*Tsolid_prev[i];
        xc[i] = z(i);
    }
}

public:
    vector_fp m_divq;
    vector_fp m_divDiffusiveSolidHeatFlux;
    vector_fp m_rad_radial;

    vector_fp m_porosity;
    vector_fp m_solidDensity;
    vector_fp m_effectiveSolidHeatConductivity;
    vector_fp m_tortuosity_factor;
    vector_fp m_solidHeatCapacity;
    vector_fp m_hydraulicDiameter;
    vector_fp m_specificArea;
    vector_fp m_extinctionCoefficient;
    vector_fp m_scatteringAlbedo;
    vector_fp m_emissivity;
    vector_fp m_hv;
    vector_fp m_insulationTransmissivity;
    vector_fp m_chemistryFactor;

    vector_fp Tgprev, Tsprev;

    std::vector<bool> m_do_solid_temperature;
    FoamStack& foams;
    double m_T_amb;
    double m_Tsolid_inlet;
    bool m_doSolidAxialRadiation;
    bool m_doSolidRadialRadiation;
    bool m_doSolid;
    bool m_solidPseudoTimeStepping;
    bool m_TsolidZeroGradient;

    bool fixedTemperature;
    void set_fixedTemperature(bool fT) { fixedTemperature = fT; }

    bool m_externTsolidEveryIteration;
    double m_relax;
    void set_firstSolve(bool b){m_firstSolve=b;}
    bool m_firstSolve;
    vector_fp xc;
    vector_fp TsFixed;
    std::string initialGuessFile;
    void setInitialGuessFile(const std::string& fname) {initialGuessFile=fname;}

    vector_fp prevMesh;
    Eigen::SparseMatrix<double,Eigen::RowMajor> mat;
    Eigen::VectorXd q;
    Eigen::VectorXd radiation_rhs;
    Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int> > solver;


    Eigen::SparseMatrix<double,Eigen::RowMajor> Tsolid_mat;
    Eigen::VectorXd Tsolid_q;
    Eigen::VectorXd Tsolid_rhs;
    Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int> > Tsolid_solver;

    void read_TsProfile()
    {
        std::ifstream read(initialGuessFile);
        std::string line;
        std::getline(read, line); //read header line
        while(std::getline(read, line))
        {
            std::stringstream Str(line);
            double xx;
            Str >> xx;
            xc.push_back(xx);

            double u;
            Str >> u;

            double t;
            Str >> t;

            double ts;
            Str >> ts;
            TsFixed.push_back(ts);
        }
    }
};

}

#endif

