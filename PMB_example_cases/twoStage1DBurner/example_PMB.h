#ifndef EXAMPLE_PMB_H
#define EXAMPLE_PMB_H

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

bool solveFlames(Cantera::Solution& sol, Cantera::FoamStack& foam_stack, const Cantera::vector_fp& mdots, const Cantera::vector_fp& phis,
        const std::string& fuel, const std::string& ox, double T0, double p0, double Tamb,
        double initialFlameProfileThickness, double totalLength, std::size_t nInitialPoints, double expected_flame_location,
        bool freeflame, bool radiation, int iterations, double residual, double relax, bool deactivateSolid, int loglevel,
        const Cantera::vector_fp& initial_Ts_x = {}, const Cantera::vector_fp& initial_Ts = {})
{
    using namespace Cantera;

    auto& gas = *sol.thermo();
    size_t nsp = gas.nSpecies();

    double pressure = p0;
    double mdot0 = mdots[0];
    double phi = phis[0];

    gas.setEquivalenceRatio(phi, fuel, ox);
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

    gas.setEquivalenceRatio(phi, fuel, ox);
    gas.setState_TP(T0, pressure);

    Cantera::PorousFlow flow(foam_stack, &gas);

    Cantera::vector_fp z(nInitialPoints+2);
    z[0]=0;

    auto dz = initialFlameProfileThickness / (nInitialPoints - 1);
    for (std::size_t iz = 0; iz < nInitialPoints; iz++)
        z[iz + 1] = iz*dz + expected_flame_location - initialFlameProfileThickness*0.5;

    z[nInitialPoints+1] = totalLength;

    flow.setupGrid(z.size(), &z[0]);

    // specify the objects to use to compute kinetic rates and transport properties

    std::unique_ptr<Cantera::Transport> trmix(Cantera::newTransportMgr("Mix", sol.thermo().get()));

    flow.setTransport(*trmix);
    flow.setKinetics(*sol.kinetics());
    flow.setPressure(pressure);
    flow.set_Tamb(Tamb);
    flow.set_Tsolid_inlet(T0);
    //flow.setTsolidZeroGradient(!radiation);

    flow.setSolid(!deactivateSolid);
    flow.setSolidAxialRadiation(radiation);
    flow.setSolidRadialRadiation(radiation);

    if (freeflame) // for freely propagating flame
    {
        flow.setFreeFlow();
    }
    else // for matrix-stabilized flame
    {
        flow.setAxisymmetricFlow();
    }


    //------- step 2: create the inlet  -----------------------

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
    //outlet.setTsolidZeroGradient(!radiation);

    //=================== create the container and insert the domains =====

    std::vector<Cantera::Domain1D*> domains { &inlet, &flow, &outlet };
    Cantera::Sim1D flame(domains);

    flame.setGridMin(1,1e-6);



    //----------- Supply initial guess----------------------

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

    flame.setMaxTimeStepCount(300);

    std::string sphi = fmt::sprintf("%1.2f",phi);
    std::string smdot= fmt::sprintf("%1.1f",mdot0);

    if (freeflame)
        flame.setFixedTemperature(0.5 * (T0 + Tad));

    if (initial_Ts_x.size() == 0)
    {
        flow.setInitialGuessFile("./data/initial_Ts_profile.txt");
        flow.read_TsProfile();
    }
    else
        flow.set_initial_guess_Ts(initial_Ts_x, initial_Ts);

    flow.set_relaxation(relax);
    //flow.setSolid(false);
    flow.set_fixedTemperature(true);
    flow.solveEnergyEqn();
    try
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

        std::vector<double> Tsprev, Tgprev;
        Cantera::vector_fp Xp(gas.nSpecies());
        bool first= true;
        for (auto tphi:phis)
        {
            for (auto mdot:mdots)
            {
                gas.setEquivalenceRatio(tphi, fuel, ox);
                gas.setState_TP(T0, pressure);
                Cantera::vector_fp x(nsp);
                gas.getMoleFractions(x.data());
                inlet.setMoleFractions(x.data());
                inlet.setMdot(mdot);


                if (first)
                {
                    flame.solve(loglevel,true);

                }
                else
                {
                    flow.setTransientTolerances(1e-2, 1e-6);
                    flow.setSteadyTolerances(1e-2, 1e-6);
                    flame.setRefineCriteria(flowdomain, ratio, slope, curve, prune);

                }

                first = false;

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
                    std::cout<<"  L2 Tg: "<<L2gas;
                    interpolate_T(gasx_curr, gasx_prev, gasTs_prev, T_tmp);
                    double L2solid = L2(gasTs_curr, T_tmp);
                    std::cout<<"  L2 Ts: "<<L2solid;

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

                std::string name = "result_phi_"+fmt::sprintf("%1.2f",tphi)+"_fuel_"+fuel+"_mdot_"+fmt::sprintf("%1.2f",mdot)+".txt";
                std::cout<<name<<std::endl;
                std::ofstream write(name);

                write << "#x Tg Ts u porosity solidDensity effectiveSolidHeatConductivity tortuosity_factor solidHeatCapacity specificArea extinctionCoefficient scatteringAlbedo emissivity hv insulationTransmissivity divq divSolidHeatConduction ";
                for (std::size_t i=0; i!=gas.nSpecies(); ++i)
                    write << gas.speciesName(i)<<" ";
                write << '\n';
                for (std::size_t i=0; i!=flow.nPoints(); ++i)
                {
                    write   << flow.grid(i)<<" "<<  flame.value(flowdomain,flow.componentIndex("T"),i) << " "<< flame.value(flowdomain,flow.componentIndex("Tsolid"),i) << " " << flame.value(flowdomain,flow.componentIndex("velocity"),i) << " "
                            << flow.m_porosity[i] << " " << flow.m_solidDensity[i] << " " << flow.m_effectiveSolidHeatConductivity[i] << " " << flow.m_tortuosity_factor[i] << " " << flow.m_solidHeatCapacity[i] << " "
                            << flow.m_specificArea[i] << " " << flow.m_extinctionCoefficient[i] << " " << flow.m_scatteringAlbedo[i] << " " << flow.m_emissivity[i] << " " << flow.m_hv[i] << " " << flow.m_insulationTransmissivity[i] << " " << flow.m_divq[i] << " " << flow.m_divDiffusiveSolidHeatFlux[i] << " ";

                    for (std::size_t k=0; k!=gas.nSpecies(); ++k)
                    {
                        write << flame.value(flowdomain,flow.componentIndex(gas.speciesName(k)),i) << " ";
                    }
                    write << '\n';
                }
                if (freeflame)
                {
                    std::cout<<"Freeflame intersitional velocity = "<<flame.value(flowdomain,flow.componentIndex("velocity"),0)<<std::endl;
                }
            }
        }
    }
    catch(const Cantera::CanteraError& e)
    {
        std::cout<<e.what()<<std::endl;
        std::cout<<e.getMessage()<<std::endl;
        std::cout<<"Info: did not converge"<<std::endl;
        return false;
    }
    catch(...)
    {
        std::cout<<"ERROR: non-Cantera exception!"<<std::endl;
        return false;
    }

    return true;

}



#endif

