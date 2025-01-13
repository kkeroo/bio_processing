import simulator, grn
import numpy as np
import matplotlib.pyplot as plt

# to test the LFSR implementation in grenmlin
def test(state):
    y2 = state[0]
    y3 = state[1]
    y4 = state[2]
    y1 = state[2] ^ state[3]

    state = [y1, y2, y3, y4]
    return state

my_grn = grn.grn()

# input species
my_grn.add_input_species("X1")
my_grn.add_input_species("X2")
my_grn.add_input_species("X3")
my_grn.add_input_species("X4")

# other species
my_grn.add_species("Y1", 0.1)
my_grn.add_species("Y2", 0.1)
my_grn.add_species("Y3", 0.1)
my_grn.add_species("Y4", 0.1)


# regulators of a gene
regulators = [{'name': 'X1', 'type': 1, 'Kd': 6, 'n': 2}]
# products of a gene
products = [{'name': 'Y2'}]

# adding a gene to the network - the first parameter specifies the rate of gene expression
my_grn.add_gene(10, regulators, products)


# regulators of a gene
regulators = [{'name': 'X2', 'type': 1, 'Kd': 5, 'n': 2}]
# products of a gene
products = [{'name': 'Y3'}]

# adding a gene to the network - the first parameter specifies the rate of gene expression
my_grn.add_gene(10, regulators, products)

# regulators of a gene
regulators = [{'name': 'X3', 'type': 1, 'Kd': 5, 'n': 2}]
# products of a gene
products = [{'name': 'Y4'}]

# adding a gene to the network - the first parameter specifies the rate of gene expression
my_grn.add_gene(15, regulators, products)


# XOR

# regulators of a gene
regulators = [{'name': 'X3', 'type': -1, 'Kd': 5, 'n': 2},
                {'name': 'X4', 'type': 1, 'Kd': 5, 'n': 3}]
# products of a gene
products = [{'name': 'Y1'}]

# adding a gene to the network - the first parameter specifies the rate of gene expression
my_grn.add_gene(10, regulators, products)

regulators = [{'name': 'X3', 'type': 1, 'Kd': 5, 'n': 2},
                {'name': 'X4', 'type': -1, 'Kd': 5, 'n': 3}]
products = [{'name': 'Y1'}]

my_grn.add_gene(10, regulators, products)

print(my_grn.genes)
print(len(my_grn.genes))

my_grn.plot_network()

def run_simulation(grn, state, n_iterations=100, filter_species=None, plot_full=False, save_plots=True):
    n_INS = len(my_grn.input_species_names)
    n_RS = len(my_grn.species_names) - n_INS
    R_0 = np.zeros(n_RS)

    full_T = None
    full_Y = None

    for i in range(n_iterations):
        T, Y = simulator.simulate_single(my_grn, state, t_end=250, R0=R_0, plot_on=False)
        if full_T is None:
            full_T = len(T)
            full_Y = Y
        else:
            full_T += len(T)
            full_Y = np.concatenate([full_Y, Y])

        last_values = Y[-1]
        X1, X2, X3, X4, Y1, Y2, Y3, Y4 = last_values

        print(f"Iteration {i}: X1: {X1}, X2: {X2}, X3: {X3}, X4: {X4}")

        R_0 = last_values[-n_RS:]

        t1 = 1 if state[0] > 50 else 0
        t2 = 1 if state[1] > 50 else 0
        t3 = 1 if state[2] > 50 else 0
        t4 = 1 if state[3] > 50 else 0
        gt_state = test([t1, t2, t3, t4])
        state = [Y1, Y2, Y3, Y4]

        t1 = 1 if Y1 > 50 else 0
        t2 = 1 if Y2 > 50 else 0
        t3 = 1 if Y3 > 50 else 0
        t4 = 1 if Y4 > 50 else 0

        testing_state = [t1, t2, t3, t4]
        
        # check if elements in testing_state and gt_state are equal
        if not testing_state == gt_state:
            print(f"Error in iteration {i}: {testing_state} != {gt_state}")
            print(f"X1: {X1}, X2: {X2}, X3: {X3}, X4: {X4}")
            print(f"Y1: {Y1}, Y2: {Y2}, Y3: {Y3}, Y4: {Y4}")
            exit(0)

    T_plot = np.arange(0, full_T)
    Y_plot = full_Y[:, filter_species] if filter_species else full_Y

    if plot_full:
        clock_period = 251

        plt.figure()
        plt.plot(T_plot, Y_plot)
        for i in range(0, len(T_plot), clock_period * 2):
            plt.axvspan(i, i + clock_period, facecolor='lightgray', alpha=0.5)

        plt.legend([my_grn.species_names[i] for i in filter_species] + ['Clock'])

        plt.xlabel('time [a.u.]')
        plt.ylabel('concentrations [a.u.]')

        plt.show()

    # Instead of ploting all four species on the same plot, we can plot them separately
    # Make a plot for each species in 4 subplots
    clock_period = 251

    for i in range(len(filter_species)):
        plt.figure()
        plt.plot(T_plot, Y_plot[:,i])
        for j in range(0, len(T_plot), clock_period * 2):
            plt.axvspan(j, j + clock_period, facecolor='lightgray', alpha=0.5)
        plt.xlabel('time [a.u.]')
        plt.ylabel(f'{my_grn.species_names[filter_species[i]]} concentration [a.u.]')
        plt.title(f'{my_grn.species_names[filter_species[i]]} concentration over time')
        plt.legend([my_grn.species_names[filter_species[i]]] + ['Clock'])

        if save_plots:
            plt.savefig(f'{my_grn.species_names[filter_species[i]]}_concentration.png')

    plt.show()


run_simulation(grn=my_grn, state=[100, 0, 0, 0], n_iterations=100, filter_species=[0, 1, 2, 3], save_plots=False)