import os

import numpy as np
import pandas as pd
import seaborn as sns

def collect_results_sweep_1(rho, theta, genome_size, sample_size, seed):
    # utilise numpy for combinations
    # https://www.kite.com/python/answers/how-to-get-all-element-combinations-of-two-numpy-arrays-in-python

    mesh_grid = np.array(np.meshgrid(rho, theta, genome_size, sample_size, seed))
    # reshape.(-1,2), -1 is unknown dimension, let numpy figure it out
    # z.reshape(-1, -1)
    # ValueError: can only specify one unknown dimension
    sweep_1_combinations = mesh_grid.T.reshape(-1, 5)

    # Load data into dataframe
    recom_est_results_dir = f"{os.getcwd()}/Output/"

    col_names = ["rho_sim", "theta_sim", "genome_size_sim", "sample_size_sim", "seed_sim", "max_rho", "max_lk"]

    df_recom_est_resutls = pd.DataFrame(columns=col_names)

    for rho, theta, genome_size, sample_size, seed in sweep_1_combinations:
        prepended_filename = f"rho_{rho}_theta_{theta}_genome_size_{int(genome_size)}_sample_size_{int(sample_size)}_seed_{int(seed)}_"
        with open(f"{recom_est_results_dir}{prepended_filename}msp_out_processed_results.csv", 'r') as results:
            results_unfiltered = results.readlines()[1].strip().split(',')

        to_append = [rho, theta, genome_size, sample_size, seed] + results_unfiltered
        df_recom_est_resutls.loc[len(df_recom_est_resutls)] = to_append

    return df_recom_est_resutls


if __name__ == '__main__':
    # Sweep 1: Recombination rate estimation
    rho_sweep_1 = [0.005, 0.0125, 0.025, 0.0375, 0.05] # unscaled r values. rho = 2 . p . N_e . r . tractlen
    theta_sweep_1 = [0.005]  # unscaled u values. theta = 2 . p . N_e . u
    genome_size_sweep_1 = [25000]
    sample_size_sweep_1 = [5,10,20,30,40,50,60,70,80,90,100]
    seed_sweep_1 = [1,2,3,4,5,6,7,8,9,10]

    recom_tract_len = 500

    collected_results_sweep_1_df = collect_results_sweep_1(rho_sweep_1, theta_sweep_1, genome_size_sweep_1,
                                                           sample_size_sweep_1, seed_sweep_1)

    # process and export df for plotting
    # Since fastsimbac does 2n*r only scaling by tract len is needed

    collected_results_sweep_1_df["scaled_rho_sim"] = collected_results_sweep_1_df["rho_sim"].apply(
        lambda x: x * recom_tract_len * 2)

    # cols_to_drop = ["rho", "sample_size", "genome_size"]
    # collected_results_sweep_1_df.drop(columns=cols_to_drop, inplace=True)

    reorder_cols = ['rho_sim', 'scaled_rho_sim', 'theta_sim', 'genome_size_sim', 'sample_size_sim',
                    'seed_sim', 'max_rho', 'max_lk']

    collected_results_sweep_1_df = collected_results_sweep_1_df.reindex(columns=reorder_cols)

    collected_results_sweep_1_df.to_csv("collected_results_msp_rho_sweep.csv", index=None)

    collected_results_sweep_1_df = collected_results_sweep_1_df.astype('float64')
    
    # Plot results
    ax = sns.boxplot(data=collected_results_sweep_1_df,x="scaled_rho_sim", y="max_rho", hue="sample_size_sim", palette="Blues")
    ax.legend(title='Genomes')
    
    ax.set_title("Rhometa_full_genome Simulated")

    ax.set(ylabel="Estimated \u03C1", xlabel="Simulated \u03C1")

    ax.figure.savefig("collected_results_msp_rho_sweep.png", dpi=500)