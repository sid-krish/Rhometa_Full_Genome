#!/usr/bin/env python
import pandas as pd
import csv

p_ij_grid = "likelihoods_grid.csv"

n_resamples = 20

df = pd.read_csv(p_ij_grid)
n_rows_df = len(df)

for i in range(n_resamples):
    # the resampled df will have the same num of rows as the original
    boot_strapped_df = df.sample(n=n_rows_df, replace=True, axis=0)
    boot_strapped_df_sums = boot_strapped_df.sum(axis=0)

    with open(f"collected_likelihoods_bootstrap_{i}.csv", 'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["rho_for_estimator", "total_of_log_likelihood"])

        for idx, val in boot_strapped_df_sums.iteritems():
            csv_writer.writerow([idx, val])

    csvfile.close()

# final results
for i in range(n_resamples):
    cl_df = pd.read_csv(f"collected_likelihoods_bootstrap_{i}.csv")

    # If multiple total_of_log_likelihoods are same, sort such that the ones connected to smallest rho_for_estimator
    # are first
    cl_df.sort_values(by=["total_of_log_likelihood", "rho_for_estimator"], ascending=[False, True], inplace=True)

    # df.reset_index(inplace=True, drop=True)

    max_rho = cl_df.iloc[0][0]
    max_lk = cl_df.iloc[0][1]

    # max_rho = round(df.iloc[0][0], 2)
    # max_lk = round(df.iloc[0][1], 2) # rounding to make it look cleaner

    cl_df.sort_values(by="rho_for_estimator", inplace=True)
    # df = df.round(2) # just so it looks cleaner
    # with open(theta_txt, 'r') as file:
    #     theta = file.readline()

    with open(f"final_results_bootstrap_{i}.txt", 'w') as file:  # open in write mode (create new file)
        file.write(f"Custom Full Genome Pairwise Recombination Rate Estimator\n")
        file.write(f"\n")

        file.write(f'-' * 56 + '\n')
        file.write(f"Max rho = {max_rho}, Max lk = {max_lk}\n")
        file.write(f'-' * 56 + '\n')

        file.write(f"\n")
        cl_df.to_csv(file, mode='a', index=None)  # open in append mode (add to new file)

    file.close()

