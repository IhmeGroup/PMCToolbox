# Code written by Thorsten Zirwes 2023
# Karlsruhe Institute of Technology
# University of Stuttgart

import cantera as ct
import numpy as np
import sys
from scipy import optimize

Tamb = 300.            # ambient temperature (used for radiation modeling) (K)
inch = 0.0254          # conversion of inches to meters
TsInit = 300.          # initial temperature of the solid (K)
constantP = ct.one_atm # gas pressure inside the porous medium (Pa)
mech = "./Stagni.yaml"   # reaction mechanism file
transmissivity = 0.6

def run_network(twoStages, solidInSecondStage, injectionInSecondStage, phi_p, fuel, mdot_p, TgasInit, Tsecond, phi_s, second_stage_pure_air, phi_g, mdot_g):

    if not twoStages:
        injectionInSecondStage = False


    if injectionInSecondStage:
        def is_none(x):
            if x is None:
                return 1
            return 0

        ecount = is_none(phi_g) + is_none(mdot_g) + is_none(mdot_p)
        if ecount != 1:
            print("ERROR! Only one of {phi_g,mdot_g} can be set!")
            sys.exit()


        if second_stage_pure_air and not phi_s is None:
            print("ERROR! Only one of {phi_s, second_stage_pure_air} can be set!")
            sys.exit()

        gas_del = ct.Solution(mech)
        gas_del.TP = TgasInit, constantP
        gas_del.set_equivalence_ratio(phi_p,fuel=fuel,oxidizer="O2:0.21,N2:0.79")
        Y_first = gas_del.Y
        if second_stage_pure_air:
            gas_del.X = "O2:0.21,N2:0.79"
        else:
            gas_del.set_equivalence_ratio(phi_s,"H2:1","O2:0.21,N2:0.79")
        Y_second = gas_del.Y

        if phi_g is None:
            Y_total = mdot_p * Y_first + (mdot_g-mdot_p)*Y_second
            gas_del.Y = Y_total
            phi_g = gas_del.equivalence_ratio()
            print("phi_g:",phi_g)

        if mdot_g is None:

            def get_phi_g(mdot_g):
                Y_total = mdot_p * Y_first + (mdot_g-mdot_p)*Y_second
                gas_del.Y = Y_total
                return gas_del.equivalence_ratio() - phi_g

            mdot_g = optimize.bisect(get_phi_g, mdot_p, 10,xtol=1e-5,rtol=1e-5)

            print("mdot_g:",mdot_g)

        if mdot_p is None:
            def get_phi_g(mdot_p):
                Y_total = mdot_p * Y_first + (mdot_g-mdot_p)*Y_second
                gas_del.Y = Y_total
                return gas_del.equivalence_ratio() - phi_g

            mdot_p = optimize.bisect(get_phi_g, 0, mdot_g,xtol=1e-5,rtol=1e-5)

            print("mdot_p:",mdot_p)

    # initial state of the gas phase containing fuel and oxidizer at the inlet
    def getInitialGas(gas,phi,fuel,T,P):
        gas.set_equivalence_ratio(phi, fuel=fuel, oxidizer="O2:0.21,N2:0.79")
        gas.TP = T,P

    # create gas phase objects for the fresh gas and burnt gas states
    gas_init = ct.Solution(mech)
    getInitialGas(gas_init,phi_p,fuel,TgasInit,constantP)

    gas_hot = ct.Solution(mech)
    getInitialGas(gas_hot,phi_p,fuel,TgasInit,constantP)
    gas_hot.equilibrate("HP")

    Y_in_inlet = gas_init.Y                # mass fractions of the inlet flow
    h_in_inlet = gas_init.enthalpy_mass    # enthalpy of the inlet flow

    h_in_second = 0.
    if injectionInSecondStage:

        gas_second = ct.Solution(mech)
        if second_stage_pure_air:
            gas_second.TPX = Tsecond, constantP, "O2:0.21,N2:0.79"
        else:
            gas_second.set_equivalence_ratio(phi_s, fuel="H2:1",oxidizer="O2:0.21,N2:0.79")

        h_in_second = gas_second.enthalpy_mass
        Y_in_second = gas_second.Y

    # Class that represents a reactor of the reactor cascade
    class ReactorProperties:
        def __init__(self, diameter, length, midpoint, chemistry, TsInit, solid, hasSolid):

            self.diameter = diameter   # diameter of the reactor (m)
            self.length = length       # length of the section (m)
            self.midpoint = midpoint   # coordinate of the reactor center (m)
            self.chemistry = chemistry # chemistry on or off in the reactor
            self.TsInit = TsInit       # initial temperature of the solid (K)
            self.solid = solid         # properties of the solid porous medium
            self.hasSolid = hasSolid   # specify whether a solid is present in the reactor
            self.mdot_total = None
            self.mdot_injection = -100.

    # class to represent the properties of the solid porous media
    class SolidProperties:
        def __init__(
                self,
                porosity,
                extinction_coefficient,
                heat_conductivity,
                specific_area,
                solid_phase,
                emissivity):

            self.porosity = porosity                   # porosity of the porous medium
            self.extinction_coefficient = extinction_coefficient         # characteristic pore diameter in m
            self.heat_conductivity = heat_conductivity # solid heat conductivity W/(m/K)
            self.specific_area = specific_area         # specific area (1/m)
            self.emissivity = emissivity               # emissivity of the solid
            self.solid_phase = solid_phase             # solution object for computing
                                                       # solid heat capacities and densities

    # Thermal conductivity of the porous medium as function of solid temperature (W/m/K)
    # Temperature-dependent fits from measurements:
    # SiC: Thermal conductivity in hot-pressed silicon carbide, D.-K. Liu, B.-W. Lin,
    #      Ceramics International, 22(5), pp. 407-414 (1996)
    def effectiveConductivitySiC(Ts):  # for silicon carbide
        return (1 - 0.86) * 1857.0 * Ts**(-0.5332)

    # YZA: Thermal conductivity of zirconia–alumina composites, N.P. Bansal, D. Zhu, 
    #      Ceramics International, 31(7), pp 911-916 (2015)
    def effectiveConductivityYZA(Ts):  # for yittria-stabilized zirconia alumina
        return 0.33

    # Define objects for the three sections of the porous media burner:
    #  YZA is the Yttria-stabilized Zirconia-Alumina section with pore density 40 PPI
    #  SiC3PPI  is a silicon carbide porous medium with a pore density of 3 PPI
    #  SiC10PPI is a silicon carbide porous medium with a pore density of 10 PPI

    YZA40PPI = SolidProperties(porosity=0.82,
                               extinction_coefficient=1340,
                               heat_conductivity=effectiveConductivityYZA,
                               specific_area=1592,
                               solid_phase=ct.Solution("YZA.yaml"),
                               emissivity=0.85)

    SiC3PPI = SolidProperties(porosity=0.86,
                               extinction_coefficient=526,
                               heat_conductivity=effectiveConductivitySiC,
                               specific_area=934,
                               solid_phase=ct.Solution("SiC.yaml"),
                               emissivity=0.85)

    SiC10PPI = SolidProperties(porosity=0.86,
                               extinction_coefficient=1340,
                               heat_conductivity=effectiveConductivitySiC,
                               specific_area=1592.,
                               solid_phase=ct.Solution("SiC.yaml"),
                               emissivity=0.85)

    SiC20PPI = SolidProperties(porosity=0.86,
                               extinction_coefficient=860,
                               heat_conductivity=effectiveConductivitySiC,
                               specific_area=1480,
                               solid_phase=ct.Solution("SiC.yaml"),
                               emissivity=0.85)

    # heat flux between the solid-phase of two neighboring reactors
    def heat_flux(left, right):
        # effective heat conductivity (lambda) in the solid

        lambda_solid = 0.
        lambda_rad = 0.
        if left.hasSolid and right.hasSolid:
            lambda_solid = 0.5 * (right.solidHeatConductivity() +
                                  left.solidHeatConductivity())

            # effective radiative heat conductivity (lambda) from the Rosseland model
            lambda_rad = 0.5 * (right.radiationConductivity() +
                                left.radiationConductivity())

        # approximate the temperature gradient between reactors
        gradT = (right.Ts - left.Ts) / (right.midpoint - left.midpoint)
        return (lambda_solid + lambda_rad) * gradT


    # New reactor class that implements equations for reacting flows in porous media (PM)
    class PMReactor(ct.ExtensibleIdealGasConstPressureReactor):
        def __init__(self, *args, props, **kwargs):
            super().__init__(*args, **kwargs)

            self.Ts = props.TsInit                     # initial temperature of the solid
            self.neighbor_left = None                  # reference to reactor on the left
            self.neighbor_right = None                 # reference to reactor on the right
            self.chemistry = props.chemistry           # turn chemistry on or off
            self.midpoint = props.midpoint             # coordinate of reactor center (m)
            self.solid = props.solid                   # solid density and heat capacity
            self.diameter = props.diameter             # reactor diameter (m)
            self.hasSolid = props.hasSolid             # is there a solid material present in the reactor?
            self.A = np.pi * (0.5 * props.diameter)**2 # cross sectional area (m^2)
            self.V = self.A * props.length             # reactor volume (m^3)
            self.molecular_weights = None
            self.species_offset = None
            self.length = props.length
            self.mdot_total = props.mdot_total
            self.mdot_injection = props.mdot_injection

        def after_initialize(self, t0):
            self.n_vars += 1                # additional equation for the solid temperature
            self.index_Ts = self.n_vars - 1
            self.species_offset = self.component_index(self.thermo.species_name(0))
            self.molecular_weights = self.thermo.molecular_weights

        def after_get_state(self, y):
            y[self.index_Ts] = self.Ts
            y[self.component_index("mass")] = self.thermo.density_mass * self.V

        # compute the effective conductivity for heat conduction within the solid
        def solidHeatConductivity(self):
            return self.solid.heat_conductivity(self.Ts)

        # effective conductivity for radiative heat transfer from the Rosseland model
        def radiationConductivity(self):
            return 16.0 * ct.stefan_boltzmann * self.Ts**3 / (3.0 * self.solid.extinction_coefficient)

        def after_update_state(self, y):
            self.Ts = y[self.index_Ts]

        # heat transfer coefficient (htc) based on a Nusselt (Nu) correlation
        def htc(self):
            d_h = 4.0 * self.solid.porosity / self.solid.specific_area # hydraulic diameter
            density = self.thermo.density_mass
            u = self.mdot_total / (density * self.solid.porosity)                   # gas velocity
            Re = density * u * d_h / self.thermo.viscosity               # Reynolds number
            lambda_gas = self.thermo.thermal_conductivity
            Pr = self.thermo.cp_mass * self.thermo.viscosity /lambda_gas # Prandtl number
            Nu = 3.7 * (Re**0.38) * (Pr**0.25)
            return Nu * self.solid.specific_area * self.thermo.thermal_conductivity / d_h

        # implement the new governing equations
        def replace_eval(self, t, LHS, RHS):

            for i in range(self.n_vars):
                LHS[i] = 1.0
                RHS[i] = 0.0

            # update solid properties with current values of the solid temperature
            self.solid.solid_phase.TP = self.Ts, constantP

            # create some variables for convenience
            Y = self.thermo.Y
            density = self.thermo.density_mass
            cp = self.thermo.cp_mass
            porosity = self.solid.porosity
            solid_density = self.solid.solid_phase.density
            solid_cp = self.solid.solid_phase.cp

            # If this reactor has no left neighbor, it is the first reactor in the cascade.
            # In this case, set initial mass fractions and enthalpies as inlet conditions
            if self.neighbor_left is None:
                Y_in = Y_in_inlet
                h_in = h_in_inlet
            else:  # otherwise, take the mass fraction and enthalpy from the left neighbor
                Y_in = self.neighbor_left.thermo.Y
                h_in = self.neighbor_left.thermo.enthalpy_mass

            # ============================================================================#
            #                          Temperature of the solid                           #
            # ============================================================================#

            if self.hasSolid:

                LHS[self.index_Ts] = (1.0 - porosity) * solid_cp * solid_density * self.V

                # heat transfer between gas and solid
                RHS[self.index_Ts] += self.htc() * (self.thermo.T - self.Ts) * self.V

                # compute the axial heat loss due to radiation to the ambience
                # if this is the first or last reactor in the cascade
                if self.neighbor_left is None or (self.neighbor_right is None and self.hasSolid) or (self.hasSolid and not self.neighbor_right.hasSolid):
                    RHS[self.index_Ts] -= ct.stefan_boltzmann * self.solid.emissivity \
                        * (self.Ts**4 - Tamb**4) * self.A

                # compute heat flux in the solid between the reactors in the cascade
                if self.neighbor_left is not None:
                    RHS[self.index_Ts] -= heat_flux(self.neighbor_left, self) * self.A
                if self.neighbor_right is not None:
                    RHS[self.index_Ts] += heat_flux(self, self.neighbor_right) * self.A

                # compute the radial radiative heat loss
                radial_radiation = ct.stefan_boltzmann * self.solid.emissivity \
                       * (self.Ts**4 - Tamb**4) * np.pi \
                       * (transmissivity*self.diameter*self.length)

                RHS[self.index_Ts] -= radial_radiation
            else:
                porosity = 1.0

            # ============================================================================#
            #                             species mass fractions                          #
            # ============================================================================#

            wdot = np.zeros(self.thermo.n_species)
            if self.chemistry:
                wdot = self.kinetics.net_production_rates # chemical source terms

            # right hand side of the mass fraction equations
            for k in range(self.thermo.n_species):
                index = k + self.species_offset
                # chemical source term
                RHS[index] += wdot[k] * self.molecular_weights[k] / density
                # convective contribution

                if injectionInSecondStage and self.mdot_injection > 0.:
                    RHS[index] += self.A / (self.V * porosity * density) * ((self.mdot_total-self.mdot_injection)*Y_in[k] + self.mdot_injection*Y_in_second[k] - self.mdot_total*Y[k])
                else:
                    RHS[index] += self.A * self.mdot_total / (self.V * porosity * density) * (Y_in[k] - Y[k])

            # ============================================================================#
            #                             Temperature of the gas                          #
            # ============================================================================#

            # right hand side of the gas phase temperature equation
            enthalpy_loss = np.dot(self.thermo.partial_molar_enthalpies
                                   / self.molecular_weights, Y_in)

            hrr = -np.dot(self.thermo.partial_molar_enthalpies, wdot)  # (W/m^3)

            Tindex = self.component_index("temperature")
            LHS[Tindex] = porosity * density * cp * self.V
            # heat transfer between solid and gas-phase
            if self.hasSolid:
                RHS[Tindex] -= self.htc() * (self.thermo.T - self.Ts) * self.V
            # convective transport
            if injectionInSecondStage and self.mdot_injection > 0.:
                RHS[Tindex] += self.A * ( (self.mdot_total-self.mdot_injection)*h_in + self.mdot_injection*h_in_second - self.mdot_total*enthalpy_loss)
            else:
                RHS[Tindex] += self.A * self.mdot_total * (h_in - enthalpy_loss)
            # chemical contribution
            RHS[Tindex] += hrr * self.V * porosity

        def before_component_index(self, name):
            if name == "Ts":
                return self.index_Ts

        def before_component_name(self, i):
            if i == self.index_Ts:
                return "Ts"

    # list of reactors that form the cascade
    reactors = []

    # Define 7 reactors, each with a different length (m), and total length of 4 inches.
    # Increasing the number of reactors improves the "resolution" and thus the accuracy

    if not twoStages:
        lengths = [
            1.5 * inch,
            0.25 * inch,
            0.25 * inch,
            # ============
            0.125 * inch,
            0.125 * inch,
            0.75 * inch,
            # ============
            0.5 * inch,
            0.5 * inch
        ]
    else:
        lengths = [
            1.5 * inch,
            0.25 * inch,
            0.25 * inch,
            # ============
            0.125 * inch,
            0.125 * inch,
            0.75 * inch,
            # ============
            0.75 * inch,
            0.125 * inch,
            0.125 * inch,
            # ============
            0.5 * inch,
            0.5 * inch,
            # ===========
            1 * inch,
            # ===========
            1 * inch
        ]
    total_length = np.sum(lengths)

    sum_lengths = 0.
    # add each reactor one by one to the list
    for rL in lengths:

        midpoint = sum_lengths + 0.5 * rL  # midpoint of the current reactor
        sum_lengths += rL

        if midpoint < 2*inch:
            # for the left half of the burner, fill all reactors with cold inert 40 PPI YZA
            props = ReactorProperties(
                diameter=2 * inch,
                length=rL,
                midpoint=midpoint,
                chemistry=False,
                TsInit=TsInit,
                solid=YZA40PPI,
                hasSolid = True)

            reactors.append(PMReactor(gas_init, props=props))

        elif midpoint < 3*inch:
            # for 2 inches < x < 3 inches, fill the reactors with hot 3 PPI SiC
            props = ReactorProperties(
                diameter=2 * inch,
                length=rL,
                midpoint=midpoint,
                chemistry=True,
                TsInit=gas_hot.T,
                solid=SiC3PPI,
                hasSolid=True)

            reactors.append(PMReactor(gas_hot, props=props))

        elif midpoint < 4*inch:
            # x >= 3 inches, fill the reactors with hot 10 PPI SiC
            props = ReactorProperties(
                diameter=2 * inch,
                length=rL,
                midpoint=midpoint,
                chemistry=True,
                TsInit=gas_hot.T,
                solid=SiC10PPI,
                hasSolid=True)

            reactors.append(PMReactor(gas_hot, props=props))

        elif midpoint < 5*inch:
            # x >= 3 inches, fill the reactors with hot 10 PPI SiC
            props = ReactorProperties(
                diameter=2 * inch,
                length=rL,
                midpoint=midpoint,
                chemistry=True,
                TsInit=gas_hot.T,
                solid=SiC3PPI,
                hasSolid=solidInSecondStage)

            reactors.append(PMReactor(gas_hot, props=props))

        elif midpoint < 6*inch:
            # x >= 3 inches, fill the reactors with hot 10 PPI SiC
            props = ReactorProperties(
                diameter=2 * inch,
                length=rL,
                midpoint=midpoint,
                chemistry=True,
                TsInit=gas_hot.T,
                solid=SiC10PPI,
                hasSolid=solidInSecondStage)

            reactors.append(PMReactor(gas_hot, props=props))

        elif midpoint < 7*inch:
            # x >= 3 inches, fill the reactors with hot 10 PPI SiC
            props = ReactorProperties(
                diameter=2 * inch,
                length=rL,
                midpoint=midpoint,
                chemistry=True,
                TsInit=gas_hot.T,
                solid=SiC20PPI,
                hasSolid=solidInSecondStage)

            reactors.append(PMReactor(gas_hot, props=props))


    # set the neighbor relations between the reactors
    for i, r in enumerate(reactors):
        if i != 0:
            r.neighbor_left = reactors[i - 1]
        if i != len(reactors) - 1:
            r.neighbor_right = reactors[i + 1]


    closest_i = -1
    distance = 1e10

    for i,r in enumerate(reactors):
        if r.midpoint < 4*inch:
            r.mdot_total = mdot_p
        elif r.midpoint > 4*inch and injectionInSecondStage:
            r.mdot_total = mdot_g
        else:
            r.mdot_total = mdot_p

        x_inlet = r.midpoint - 0.5*r.length
        d = np.abs(x_inlet - 4*inch)
        if d < distance:
            distance = d
            closest_i = i

    if injectionInSecondStage:
        reactors[closest_i].mdot_injection = mdot_g - mdot_p

    # create the cascade and set the numerical tolerances for time integration
    net = ct.ReactorNet(reactors)
    net.max_steps = 100000
    net.atol = 1e-11
    net.rtol = 1e-6

    net.atol_sensitivity = 1e-3
    net.rtol_sensitivity = 0.05

    # simulate the flow through burner until the steady-state
    #net.advance_to_steady_state()
    net.advance(net.time+3600.)
    return reactors

def pollutants(gas):
    X = gas.X
    X[gas.species_index("H2O")] = 0.
    gas.X = X
    X = gas.X
    NOx = gas.X[gas.species_index("NO")]*5.95/(20.95-X[gas.species_index("O2")])*1e6;
    NH3 = X[gas.species_index("NH3")]*5.95/(20.95-X[gas.species_index("O2")])*1e6;
    H2 = 1e2 * X[gas.species_index("H2")];
    N2O = (X[gas.species_index("N2O")])*5.95/(20.95-X[gas.species_index("O2")])*1e6;
    return NOx, NH3, H2, N2O

import matplotlib.pyplot as plt

reactors = run_network(twoStages=True, solidInSecondStage=True, injectionInSecondStage=True, phi_p=1.3, fuel="NH3:0.7,H2:0.3",
                        mdot_p=0.3, TgasInit=300, Tsecond=300, phi_s=None, second_stage_pure_air=True, phi_g=None, mdot_g=0.4)
x = [r.midpoint for r in reactors]
Tg = [r.T for r in reactors]
Ts = [r.Ts for r in reactors]

plt.plot(x,Tg,label="Tg")
plt.plot(x,Ts,label="Ts")
plt.xlabel("x (m)")
plt.ylabel("T (K)")
plt.legend()
plt.grid(True)
plt.savefig("T.png")

