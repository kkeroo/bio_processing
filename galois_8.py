import simulator, grn
import numpy as np
import matplotlib.pyplot as plt


my_grn = grn.grn()

delta = 0.1
alpha = 10

# input species
my_grn.add_input_species("X1")
my_grn.add_input_species("X2")
my_grn.add_input_species("X3")
my_grn.add_input_species("X4")
my_grn.add_input_species("X5")
my_grn.add_input_species("X6")
my_grn.add_input_species("X7")
my_grn.add_input_species("X8")

# other species
my_grn.add_species("Y1", delta=delta)
my_grn.add_species("Y2", delta=delta)
my_grn.add_species("Y3", delta=delta)
my_grn.add_species("Y4", delta=delta)
my_grn.add_species("Y5", delta=delta)
my_grn.add_species("Y6", delta=delta)
my_grn.add_species("Y7", delta=delta)
my_grn.add_species("Y8", delta=delta)


# Always zero
regulators = [{'name': 'X8', 'type': 1, 'Kd': 1, 'n': 2}]
# products of a gene
products = [{'name': 'Y1'}]

# adding a gene to the network - the first parameter specifies the rate of gene expression
my_grn.add_gene(alpha, regulators, products)

# regulators of a gene
regulators = [{'name': 'X1', 'type': 1, 'Kd': 5, 'n': 2}]
# products of a gene
products = [{'name': 'Y2'}]

# adding a gene to the network - the first parameter specifies the rate of gene expression
my_grn.add_gene(alpha, regulators, products)

# regulators of a gene
regulators = [{'name': 'X2', 'type': 1, 'Kd': 5, 'n': 2}]
# products of a gene
products = [{'name': 'Y3'}]

# adding a gene to the network - the first parameter specifies the rate of gene expression
my_grn.add_gene(alpha, regulators, products)

# 4 = (8 ∧ ¬3) ∨ (¬8 ∧ 3)
regulators = [{'name': 'X8', 'type': 1, 'Kd': 1, 'n': 2},
                {'name': 'X3', 'type': -1, 'Kd': 1, 'n': 2}]
products = [{'name': 'Y4'}]
my_grn.add_gene(10, regulators, products)

regulators = [{'name': 'X8', 'type': -1, 'Kd': 1, 'n': 2},
                {'name': 'X3', 'type': 1, 'Kd': 1, 'n': 2}]
products = [{'name': 'Y4'}]
my_grn.add_gene(10, regulators, products)

# regulators of a gene
regulators = [{'name': 'X4', 'type': 1, 'Kd': 5, 'n': 2}]
# products of a gene
products = [{'name': 'Y5'}]

# adding a gene to the network - the first parameter specifies the rate of gene expression
my_grn.add_gene(alpha, regulators, products)

# 6 =  (8 ∧ ¬5) ∨ (¬8 ∧ 5)
regulators = [{'name': 'X8', 'type': 1, 'Kd': 1, 'n': 2},
                {'name': 'X5', 'type': -1, 'Kd': 1, 'n': 2}]
products = [{'name': 'Y6'}]
my_grn.add_gene(10, regulators, products)

regulators = [{'name': 'X8', 'type': -1, 'Kd': 1, 'n': 2},
                {'name': 'X5', 'type': 1, 'Kd': 1, 'n': 2}]
products = [{'name': 'Y6'}]
my_grn.add_gene(10, regulators, products)

# regulators of a gene
regulators = [{'name': 'X6', 'type': 1, 'Kd': 5, 'n': 2}]
# products of a gene
products = [{'name': 'Y7'}]

# adding a gene to the network - the first parameter specifies the rate of gene expression
my_grn.add_gene(alpha, regulators, products)

# 8 = (8 ∧ ¬7) ∨ (¬8 ∧ 7)
regulators = [{'name': 'X8', 'type': 1, 'Kd': 1, 'n': 2},
                {'name': 'X7', 'type': -1, 'Kd': 1, 'n': 2}]
products = [{'name': 'Y8'}]
my_grn.add_gene(10, regulators, products)

regulators = [{'name': 'X8', 'type': -1, 'Kd': 1, 'n': 2},
                {'name': 'X7', 'type': 1, 'Kd': 1, 'n': 2}]
products = [{'name': 'Y8'}]
my_grn.add_gene(10, regulators, products)

print(my_grn.genes)
print(len(my_grn.genes))

my_grn.plot_network()

def run_simulation(grn, state, n_iterations=100, filter_species=None, plot_full=False, save_plots=True, print_output=False):
    n_INS = len(grn.input_species_names)
    n_RS = len(grn.species_names) - n_INS
    R_0 = np.zeros(n_RS)

    full_T = None
    full_Y = None

    for i in range(n_iterations):
        T, Y = simulator.simulate_single(grn, state, t_end=250, R0=R_0, plot_on=False)
        if full_T is None:
            full_T = len(T)
            full_Y = Y
        else:
            full_T += len(T)
            full_Y = np.concatenate([full_Y, Y])

        last_values = Y[-1]
        X1, X2, X3, X4, X5, X6, X7, X8, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8 = last_values

        if print_output:
            print(f"Iteration {i}: X1: {X1:.1f}, X2: {X2:.1f}, X3: {X3:.1f}, X4: {X4:.1f} X5: {X5:.1f} X6: {X6:.1f} "
                  f"X7: {X7:.1f} X8: {X8:.1f}")

        R_0 = last_values[-n_RS:]

        state = [Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8]

    T_plot = np.arange(0, full_T)
    Y_plot = full_Y[:, filter_species] if filter_species else full_Y

    if plot_full:
        clock_period = 251

        plt.figure()
        plt.plot(T_plot, Y_plot)
        for i in range(0, len(T_plot), clock_period * 2):
            plt.axvspan(i, i + clock_period, facecolor='lightgray', alpha=0.5)

        plt.legend([grn.species_names[i] for i in filter_species] + ['Clock'])

        plt.xlabel('time [a.u.]')
        plt.ylabel('concentrations [a.u.]')

        plt.show()

    # Instead of ploting all four species on the same plot, we can plot them separately
    # Make a plot for each species in 4 subplots
    clock_period = 251

    for i in range(len(filter_species)):
        plt.figure()
        T_plot = np.arange(0, full_T)
        plt.plot(T_plot, full_Y[:,8+i])
        for j in range(0, len(T_plot), clock_period * 2):
            plt.axvspan(j, j + clock_period, facecolor='lightgray', alpha=0.5)
        plt.xlabel('time [a.u.]')
        plt.ylabel(f'{grn.species_names[filter_species[i]]} concentration [a.u.]')
        plt.title(f'{grn.species_names[filter_species[i]]} concentration over time')
        plt.legend([grn.species_names[filter_species[i]]] + ['Clock'])

        if save_plots:
            plt.savefig(f'{grn.species_names[filter_species[i]]}_concentration.png')

    plt.show()

run_simulation(grn=my_grn, state=[100, 0, 0, 0, 0, 0, 0, 0], n_iterations=50, filter_species=[0, 1, 2, 3, 4, 5, 6, 7], print_output=True, save_plots=False)
