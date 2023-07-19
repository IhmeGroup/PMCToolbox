//! @file StFlow_PMB.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

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


#include "cantera/oneD/StFlow.h"
#include "cantera/oneD/StFlow_PMB.h"
#include "cantera/oneD/Boundary1D_PMB.h"
#include "cantera/oneD/refine.h"
#include "cantera/base/ctml.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/numerics/funcs.h"
#include "cantera/base/global.h"
#include <iomanip>
#include <set>
#include <fstream>
using namespace std;

namespace Cantera
{

    PorousFlow::PorousFlow(FoamStack& foamList, ThermoPhase* ph, size_t nsp, size_t points)
        :   StFlow(ph, nsp, points),
        foams(foamList),
        m_doSolidAxialRadiation(true),
        m_doSolidRadialRadiation(true),
        m_doSolid(true),
        m_solidPseudoTimeStepping(false),
        m_TsolidZeroGradient(false),
        m_externTsolidEveryIteration(true),
        m_relax(1.0),
        m_firstSolve(true)
    {
        m_type = cAxisymmetricStagnationFlow;
        m_dovisc = false;
        setBounds(c_offset_Tsolid, 0, 1.0e20);
        m_refiner->setActive(c_offset_Tsolid, false);
        m_do_solid_temperature.resize(m_points,false);
    }
    void PorousFlow::solveEnergyEqn(size_t j)
    {
        StFlow::solveEnergyEqn(j);
        solve_solid_temperature(j);
    }
    void PorousFlow::fixTemperature(size_t j)
    {
        StFlow::fixTemperature(j);
        bool changed = false;
        if (j == npos) {
            for (size_t i = 0; i < m_points; i++) {
                if (m_do_solid_temperature[i]) {
                    changed = true;
                }
                m_do_solid_temperature[i] = false;
            }
        } else {
            if (m_do_solid_temperature[j]) {
                changed = true;
            }
            m_do_solid_temperature[j] = false;
        }
        m_refiner->setActive(c_offset_Tsolid, false);
        if (changed) {
            needJacUpdate();
        }
    }
    void PorousFlow::solve_solid_temperature(size_t j)
    {
        bool changed = false;
        if (j == npos) {
            for (size_t i = 0; i < m_points; i++) {
                if (!m_do_solid_temperature[i]) {
                    changed = true;
                }
                m_do_solid_temperature[i] = true;
            }
        } else {
            if (!m_do_solid_temperature[j]) {
                changed = true;
            }
            m_do_solid_temperature[j] = true;
        }
        m_refiner->setActive(c_offset_U, true);
        m_refiner->setActive(c_offset_V, true);
        m_refiner->setActive(c_offset_T, true);
        m_refiner->setActive(c_offset_Tsolid, true);
        if (changed) {
            needJacUpdate();
        }
    }
} // namespace

