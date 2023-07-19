#ifndef CT_BOUNDARY1D_PMB_H
#define CT_BOUNDARY1D_PMB_H

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

#include "Boundary1D.h"
#include "StFlow_PMB.h"

namespace Cantera
{

class InletWithSolid1D : public Inlet1D
{
    public:
        InletWithSolid1D()
            : Inlet1D(),
            m_Tsolid(-1.0),
            foams(nullptr)
    {
    }


    virtual void setFoams(FoamStack* p) {foams=p;}
    virtual void setTsolid(double Tsolid)
    {
        m_Tsolid = Tsolid;
        needJacUpdate();
    }
    virtual double getTsolid() const
    {
        return m_Tsolid;
    }

    virtual void init()
    {
        Inlet1D::init();
    }

    void eval(size_t jg, double* xg, double* rg,
            integer* diagg, double rdt)
    {
        if (jg != npos && (jg + 2 < firstPoint() || jg > lastPoint() + 2)) {
            return;
        }

        double porosity = foams->get_local_foam_properties(0.0,m_Tsolid).porosity;

        if (m_ilr == LeftInlet) {
            // Array elements corresponding to the first point of the flow domain
            double* xb = xg + m_flow->loc();
            double* rb = rg + m_flow->loc();


            if (!m_flow->fixed_mdot())
            {
                //rb[c_offset_Tsolid] -= m_Tsolid; //causes 2*Ts(0) artifact at inlet
            }


            // The first flow residual is for u. This, however, is not modified by
            // the inlet, since this is set within the flow domain from the
            // continuity equation.

            // spreading rate. The flow domain sets this to V(0),
            // so for finite spreading rate subtract m_V0.
            rb[c_offset_V] -= m_V0;

            if (m_flow->doEnergy(0)) {
                // The third flow residual is for T, where it is set to T(0).  Subtract
                // the local temperature to hold the flow T to the inlet T.
                rb[c_offset_T] -= m_temp;
            }

            if (m_flow->fixed_mdot()) {
                // The flow domain sets this to -rho*u. Add mdot to specify the mass
                // flow rate.
                rb[c_offset_L] += m_mdot;
            } else {
                // if the flow is a freely-propagating flame, mdot is not specified.
                // Set mdot equal to rho*u, and also set lambda to zero.
                m_mdot = m_flow->density(0)*xb[0];
                rb[c_offset_L] = xb[c_offset_L];
            }

            // add the convective term to the species residual equations
            for (size_t k = 0; k < m_nsp; k++) {
                if (k != m_flow_right->leftExcessSpecies()) {

                    if (m_flow->fixed_mdot()) {
                        rb[c_offset_Y+k] += m_mdot*m_yin[k];
                    }
                    else {
                        rb[c_offset_Y+k] += porosity*m_mdot*m_yin[k];
                    }
                }
            }

        } else {
            //note: currently, only the left boundary can act as a porous medium inlet

            // Array elements corresponding to the last point in the flow domain
            double* rb = rg + loc() - m_flow->nComponents();
            rb[c_offset_V] -= m_V0;
            if (m_flow->doEnergy(m_flow->nPoints() - 1)) {
                rb[c_offset_T] -= m_temp;
            }
            rb[c_offset_U] += m_mdot;
            for (size_t k = 0; k < m_nsp; k++) {
                if (k != m_flow_left->rightExcessSpecies()) {
                    rb[c_offset_Y+k] += m_mdot * m_yin[k];
                }
            }
        }
    }
    // Set the total mass flow rate.
    virtual void setMdot(double mdot) {
        m_mdot = mdot;
    }

    protected:
        double m_Tsolid;
        FoamStack* foams;
};

class OutletWithSolid1D : public Outlet1D
{
public:
    OutletWithSolid1D() : Outlet1D(), m_TsolidZeroGradient(false), Tamb(-1.0), foams(nullptr) {}

    void eval(size_t jg, double* xg, double* rg, integer* diagg,
            double rdt)
    {
        if (jg != npos && (jg + 2 < firstPoint() || jg > lastPoint() + 2)) {
            return;
        }

        // start of local part of global arrays
        double* x = xg + loc();
        double* r = rg + loc();
        integer* diag = diagg + loc();

        if (m_flow_right) {
            size_t nc = m_flow_right->nComponents();
            double* xb = x;
            double* rb = r;
            rb[c_offset_U] = xb[c_offset_L];
            if (m_flow_right->doEnergy(0)) {
                rb[c_offset_T] = xb[c_offset_T] - xb[c_offset_T + nc];
            }
            for (size_t k = c_offset_Y; k < nc; k++) {
                rb[k] = xb[k] - xb[k + nc];
            }
        }

        if (m_flow_left) {
            size_t nc = m_flow_left->nComponents();
            double* xb = x - nc;
            double* rb = r - nc;
            int* db = diag - nc;

            // zero Lambda
            if (m_flow_left->fixed_mdot()) {
                rb[c_offset_U] = xb[c_offset_L];
            }

            if (m_flow_left->doEnergy(m_flow_left->nPoints()-1)) {
                rb[c_offset_T] = xb[c_offset_T] - xb[c_offset_T - nc]; // zero T gradient
            }
            size_t kSkip = c_offset_Y + m_flow_left->rightExcessSpecies();
            for (size_t k = c_offset_Y; k < nc; k++) {
                if (k != kSkip) {
                    rb[k] = xb[k] - xb[k - nc]; // zero mass fraction gradient
                    //db[k] = 0;
                }
            }
        }

        if (m_flow_right) {
            std::cerr<<"Not implemented"<<std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (m_flow_left) {
            size_t nc = m_flow_left->nComponents();
            double* xb = x - nc;
            double* rb = r - nc;
            int* db = diag - nc;
            double L = m_flow_left->grid(m_flow_left->nPoints()-1);
            double Lm1 = m_flow_left->grid(m_flow_left->nPoints()-2);
            const auto& foam = foams->get_local_foam_properties(L, xb[c_offset_Tsolid]);
            double porosity =foam.porosity;
            double lambda=foam.effective_heat_conductivity;
            double emissivity=foam.emissivity;
            if (m_TsolidZeroGradient)
            {
                rb[c_offset_Tsolid] = xb[c_offset_Tsolid] - xb[c_offset_Tsolid - nc];
                //db[c_offset_Tsolid] = 0;
            }
            else
            {
                rb[c_offset_Tsolid] = (xb[c_offset_Tsolid] - xb[c_offset_Tsolid - nc])*lambda/(L-Lm1) + emissivity*StefanBoltz*(1.-porosity)*(std::pow(xb[c_offset_Tsolid],4.0) - std::pow(Tamb,4.0));
            }
        }
    }

    virtual void setTsolidZeroGradient(bool flag) {m_TsolidZeroGradient = flag;}
    virtual bool getTsolidZeroGradient() const {return m_TsolidZeroGradient;}
    virtual void setFoams(FoamStack* p) {foams=p;}
    virtual void setTamb(double Ta) {Tamb=Ta;}
protected:
    bool m_TsolidZeroGradient;
    double Tamb;
    FoamStack* foams;
};

}
#endif
