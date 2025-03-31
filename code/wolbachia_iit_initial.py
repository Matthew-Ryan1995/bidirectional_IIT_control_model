#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 14:09:26 2024

1. there is a bug
2. Need general parameter input
3. Need secondary parameter calculations
4. Need CI
5. Need interventions: Events in gillespy2
    a+ wPip population
7. Need Frieds Index
6. Need immigration

@author: rya200
"""

import matplotlib.pyplot as plt
import gillespy2
import numpy as np
from scipy.integrate import solve_ivp
from datetime import datetime


start_time = datetime.now()

# %%

params = {
    "numMaleClasses": 20,
    "numFemaleClasses": 1,
    "numImmatureClasses": 10,
    "carryingCapacityByPopulation_vector": 100,
    "transitionRate_imm": 1,
    "transitionRate_m": 1,
    "transitionRate_f": 1,
    "maleDeathRate": 1/7.8,
    "femaleDeathRate": 1/10,
    "propMated": 0.8,
    "femaleBirthRate": 0.253,
    "propFemale": 0.5,
    "propMale": 0.5
}


def get_secondary_params(params):
    theta = params["maleDeathRate"] * params["propFemale"] / \
        (params["femaleDeathRate"] * params["propMale"])
    phi = params["transitionRate_m"] / \
        (params["transitionRate_m"] + params["maleDeathRate"])

    F_bar = params["carryingCapacityByPopulation_vector"] * \
        theta * (1 - params["propMated"]) / (1 + theta)
    F_bar_mated = params["carryingCapacityByPopulation_vector"] * \
        theta * (params["propMated"]) / (1 + theta)
    M_bar = params["carryingCapacityByPopulation_vector"] / (1 + theta)
    M1 = M_bar / ((1 - phi**(params["numMaleClasses"] - 1)) / (1-phi) + (
        params["transitionRate_m"] * phi**(params["numMaleClasses"] - 2))/(params["maleDeathRate"]))

    I_bar = params["transitionRate_m"] * M1 / \
        (phi * params["transitionRate_m"] * params["propMale"])
    I_total = params["numImmatureClasses"] * I_bar
    I_max = I_total / (1 - (params["transitionRate_imm"]
                       * I_bar) / (params["femaleBirthRate"] * F_bar_mated))

    matingRate = (((params["transitionRate_imm"] * params["propFemale"]) *
                  I_bar - params["femaleDeathRate"] * F_bar) / (F_bar * M_bar)) * M_bar

    secondary_params = {
        "theta": theta,
        "phi": phi,
        "F_bar": F_bar,
        "F_bar_mated": F_bar_mated,
        "M_bar": M_bar,
        "M1": M1,
        "I_bar": I_bar,
        "I_total": I_total,
        "I_max": I_max,
        "matingRate": matingRate
    }

    return secondary_params


def get_steady_state(params):

    ss = {}

    for idx in range(params["numImmatureClasses"]):
        ss[f"I{idx + 1}"] = np.ceil(params["I_bar"])

    # todo: check in num male class = 1, treat different
    for idx in range(params["numMaleClasses"]):
        if (idx + 1) < params["numMaleClasses"]:
            ss[f"M{idx+1}"] = np.ceil(params["phi"]**(idx) * params["M1"])
        else:
            ss[f"M{idx+1}"] = np.ceil(params["M1"] * (params["transitionRate_m"]
                                      * params["phi"]**(idx-1))/params["maleDeathRate"])

    ss["F_unmated"] = np.ceil(params["F_bar"])

    for idx in range(params["numMaleClasses"]):
        ss[f"F_mated_{idx+1}"] = np.ceil(params["matingRate"] * params["F_bar"]
                                         * ss[f"M{idx+1}"] / (params["femaleDeathRate"] * params["M_bar"]))
    return ss


# %%
def iit(parameter_values, ss, t_end=20):

    model = gillespy2.Model(name="iit_programme")

    # Parameters

    # Number of immature and adult male compartments
    k = parameter_values["numImmatureClasses"]
    K = parameter_values["numMaleClasses"]

    # Rate of transition between compartments
    gamma = gillespy2.Parameter(
        name="gamma", expression=parameter_values["transitionRate_imm"])
    sigma = gillespy2.Parameter(
        name="sigma", expression=parameter_values["transitionRate_m"])

    mu_M = gillespy2.Parameter(
        name="mu_M", expression=parameter_values["maleDeathRate"])
    mu_F = gillespy2.Parameter(
        name="mu_F", expression=parameter_values["femaleDeathRate"])

    lam = gillespy2.Parameter(
        name="lam", expression=parameter_values["femaleBirthRate"])

    eta = gillespy2.Parameter(
        name="eta", expression=parameter_values["matingRate"])

    I_max = gillespy2.Parameter(
        name="I_max", expression=parameter_values["I_max"])

    model.add_parameter([gamma, sigma, mu_M, mu_F, lam, eta, I_max])

    # Compartments
    compartment_list = []
    I_compartments = []
    I_tot = 0
    for idx in range(k):
        I_compartments.append(gillespy2.Species(
            name=f"I{idx+1}", initial_value=ss[f"I{idx+1}"]))
        I_tot += ss[f"I{idx+1}"]
    compartment_list.append(I_compartments)

    M_compartments = []
    M_compartments.append(gillespy2.Species(name="M1", initial_value=ss["M1"]))
    M_tots = ss["M1"]
    for idx in range(1, K):
        M_compartments.append(gillespy2.Species(
            name=f"M{idx+1}", initial_value=ss[f"M{idx+1}"]))
        M_tots += ss[f"M{idx+1}"]

    M_compartments.append(gillespy2.Species(
        name="M_totals", initial_value=M_tots))

    compartment_list.append(M_compartments)

    compartment_list.append(gillespy2.Species(
        name="F_unmated", initial_value=ss["F_unmated"]))

    F_mated_compartments = []
    for idx in range(K):
        F_mated_compartments.append(gillespy2.Species(
            name=f"F_mated_{idx+1}", initial_value=ss[f"F_mated_{idx+1}"]))
    compartment_list.append(F_mated_compartments)

    compartment_list.append(gillespy2.Species(
        name="I_total", initial_value=I_tot))

    model.add_species(compartment_list)

    # Transition rates
    rates = []
    # Immature
    for idx in range(k-1):
        rates.append(gillespy2.Reaction(name=f"I{idx+1}_age",
                                        reactants={f"I{idx+1}": 1},
                                        products={f"I{idx+2}": 1},
                                        rate=gamma))
    rates.append(gillespy2.Reaction(name=f"I{k}_age_to_male",
                                    reactants={f"I{k}": 1, "I_total": 1},
                                    products={f"M1": 1, "M_totals": 1},
                                    propensity_function=f"I{k} * gamma/2"))
    rates.append(gillespy2.Reaction(name=f"I{k}_age_to_female",
                                    reactants={f"I{k}": 1, "I_total": 1},
                                    products={f"F_unmated": 1},
                                    propensity_function=f"I{k} * gamma/2"))

    # Male aging
    for idx in range(K-1):
        rates.append(gillespy2.Reaction(name=f"M{idx+1}_age",
                                        reactants={f"M{idx+1}": 1},
                                        products={f"M{idx+2}": 1},
                                        rate=sigma))

    # Male deaths
    for idx in range(K):
        rates.append(gillespy2.Reaction(name=f"M_{idx+1}_death",
                                        reactants={
                                            f"M{idx+1}": 1, "M_totals": 1},
                                        products={},
                                        propensity_function=f"M{idx+1} * mu_M"))

    # Female deaths
    rates.append(gillespy2.Reaction(name="F_unmated_death",
                                    reactants={"F_unmated": 1},
                                    products={},
                                    rate=mu_F))
    for idx in range(K):
        rates.append(gillespy2.Reaction(name=f"F_(mated, {idx+1})_death",
                                        reactants={f"F_mated_{idx+1}": 1},
                                        products={},
                                        rate=mu_F))

    # Mating
    for idx in range(K):
        rates.append(gillespy2.Reaction(name=f"F_mate_with_M{idx+1}",
                                        reactants={"F_unmated": 1},
                                        products={f"F_mated_{idx+1}": 1},
                                        propensity_function=f"eta * F_unmated * M{idx+1} / M_totals"))

    # Births
    for idx in range(K):
        rates.append(gillespy2.Reaction(name=f"F_(mated, {idx+1})_brith",
                                        reactants={},
                                        products={"I1": 1, "I_total": 1},
                                        propensity_function=f"F_mated_{idx+1} * lam * (I_max - I_total)/I_max"))

    model.add_reaction(rates)

    tspan = gillespy2.TimeSpan.linspace(t=t_end, num_points=t_end + 1)
    model.timespan(tspan)

    return model

# %%


# %%

s_params = get_secondary_params(params)

pp = params.copy()
pp.update(s_params)

ss = get_steady_state(pp)


num_trajectories = 100

M = iit(parameter_values=pp, ss=ss, t_end=1000)
results = M.run(number_of_trajectories=num_trajectories)

# %%
k = 10
K = 20

plt.figure()
for index in range(0, num_trajectories):
    trajectory = results[index]

    immatures = trajectory["I1"]
    for idx in range(1, k):
        immatures += trajectory[f"I{idx+1}"]

    plt.plot(trajectory['time'], immatures, 'g', alpha=0.1)
plt.title("Immatures")
plt.xlabel("Days")
plt.ylabel("Number")
plt.show()

plt.figure()
for index in range(0, num_trajectories):
    trajectory = results[index]

    males = trajectory["M1"]
    for idx in range(1, K):
        males += trajectory[f"M{idx+1}"]

    plt.plot(trajectory['time'], males,   'r', alpha=0.1)

plt.title("Males")
plt.xlabel("Days")
plt.ylabel("Number")
plt.show()

plt.figure()
for index in range(0, num_trajectories):
    trajectory = results[index]
    females = trajectory["F_unmated"]
    for idx in range(K):
        females += trajectory[f"F_mated_{idx+1}"]

    plt.plot(trajectory['time'], females,   'b', alpha=0.1)
plt.title("Females")
plt.xlabel("Days")
plt.ylabel("Number")
plt.show()

# %%

print(f"Time taken: {datetime.now()-start_time}")
