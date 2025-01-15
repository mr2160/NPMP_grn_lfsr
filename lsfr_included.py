import numpy as np
import simulator, grn
import matplotlib.pyplot as plt
from components import *

# -----------------------------------------------------
def add_lfsr(my_grn, name_prefix, clock_species, kd=10, n=4, delta=0.7, alpha=40):
    """
    Adds a 3-bit LFSR to 'my_grn' using D flip-flops already implemented.
    """
    
    feedback_species = f"{name_prefix}_feedback"
    my_grn.add_species(feedback_species, delta=delta)
    
    dff1 = f"{name_prefix}1"
    dff2 = f"{name_prefix}2"
    dff3 = f"{name_prefix}3"
    
    add_d_flip_flop(my_grn, dff1, feedback_species, clock_species, kd, n, delta, alpha)
    add_d_flip_flop(my_grn, dff2, f"{dff1}_DL2_Q", clock_species, kd, n, delta, alpha)
    add_d_flip_flop(my_grn, dff3, f"{dff2}_DL2_Q", clock_species, kd, n, delta, alpha)
    
    tap1 = f"{dff2}_DL2_Q"
    tap2 = f"{dff3}_DL2_Q"
    

    my_grn.add_gene(
        alpha,
        [
            {'name': tap1, 'type': -1, 'Kd': kd, 'n': n},
            {'name': tap2, 'type':  1, 'Kd': kd, 'n': n},
        ],
        [{'name': feedback_species}]
    )
    my_grn.add_gene(
        alpha,
        [
            {'name': tap1, 'type':  1, 'Kd': kd, 'n': n},
            {'name': tap2, 'type': -1, 'Kd': kd, 'n': n},
        ],
        [{'name': feedback_species}]
    )


# -----------------------------------------------------
# 2) Build an LFSR GRN and simulate:

if __name__ == "__main__":
    lfsr_grn = grn.grn()
    
    lfsr_grn.add_input_species("C")

    # LFSR parameters
    alpha = 40
    delta = 0.7
    kd = 10
    n = 4

    add_lfsr(lfsr_grn, "LFSR", "C", kd=kd, n=n, delta=delta, alpha=alpha)

    lfsr_grn.set_species_ic("LFSR1_DL2_Q", 100.0)  # for example
    lfsr_grn.set_species_ic("LFSR2_DL2_Q", 0.0)
    lfsr_grn.set_species_ic("LFSR3_DL2_Q", 100.0)


    inputs = [(100,), (0,)] * 10

    # Simulate
    T, Y = simulator.simulate_sequence(
        lfsr_grn,
        inputs,
        t_single=250,
        plot_on=False,  # We'll manually plot below
        displayed=[]
    )

    displayed = [
        'C',
        'LFSR_feedback',
        'LFSR1_DL2_Q',
        'LFSR2_DL2_Q',
        'LFSR3_DL2_Q'
    ]

    num_subplots = len(displayed)
    fig, axs = plt.subplots(num_subplots, 1, figsize=(10, 2.5 * num_subplots), sharex=True)

    for i, name in enumerate(displayed):
        idx = lfsr_grn.species_names.index(name)
        axs[i].plot(T, Y[:, idx], label=name, linewidth=2)
        axs[i].set_ylabel(name)
        axs[i].legend(loc='upper right')

    axs[-1].set_xlabel("Time")
    plt.tight_layout()
    plt.show()
