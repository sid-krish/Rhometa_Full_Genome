#! /usr/bin/env python
import sys

run_result = sys.argv[1]

with open(run_result, 'r') as f:
    lines = f.readlines()

max_rho, max_lk = lines[4].split(',')

max_rho = max_rho.strip()
max_lk = max_lk.strip()

max_rho = max_rho.strip("Max rho=")
max_lk = max_lk.strip("Max lk=")

with open("processed_results.csv", 'w') as file:  # open in write mode (create new file)
    file.write(f"max_rho,max_lk\n")
    file.write(f"{max_rho},{max_lk}\n")

