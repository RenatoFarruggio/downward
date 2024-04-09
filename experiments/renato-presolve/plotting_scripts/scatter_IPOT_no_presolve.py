import json

#print("Hello, world!")

file_path = '2024-02-06_4_no_preprocessing-eval/properties'

with open(file_path, 'r') as file:
    data = json.load(file)

data_points = []

ignored_algorithms = ["config_DEL_default",
                      "config_DEL_no_preprocessing",
                      "config_DPOT_default",
                      "config_DPOT_no_preprocessing",
                      #"config_IPOT_default",
                      #"config_IPOT_no_preprocessing",
                      "config_OCP_default",
                      "config_OCP_no_preprocessing",
                      "config_PHO_default",
                      "config_PHO_no_preprocessing",
                      "config_SEH_default",
                      "config_SEH_no_preprocessing",
                      ]

for key in data:
    #print(f"{key=}")
    error = data[key].get("error", None)
    #print(f"{error=}")
    if error is None:
        print(f"Error: No error found for {key=}")
        exit()
    if error != "success":
        error = "failed"

    id = data[key].get("id", None)
    if id is None:
        print(f"Error: No id found for {key=}")
        exit()
    domain, prob = id[1], id[2]

    #print(f"{domain=}")
    #print(f"{prob=}")

    algorithm = data[key].get("algorithm", None)
    algorithm = algorithm.split(':')[-1]
    if algorithm is None:
        print(f"Error: No algorithm found for {key=}")
        exit()
    

    if algorithm in ignored_algorithms:
        continue

    if algorithm == "config_IPOT_default":
        algorithm = "default"
    if algorithm == "config_IPOT_no_preprocessing":
        algorithm = "barrier"

    total_time = data[key].get("total_time", None)
    if total_time is None and error == "success":
        print(f"Error: No total_time found for {key=}")
        exit()

    if error != "success":
        total_time = 300

    data_points.append((algorithm, f"{domain}:{prob}", total_time, error))







import pandas as pd
import matplotlib.pyplot as plt

df = pd.DataFrame(data_points, columns=['Algorithm', 'Instance', 'Time', 'Status'])

# Pivot the DataFrame to get matching times for each instance
df_pivoted = df.pivot(index='Instance', columns='Algorithm', values='Time').reset_index()



plt.figure(figsize=(8, 6))

plt.scatter(df_pivoted['default'], df_pivoted['barrier'], c='black', edgecolors='None', marker='.')

plt.xscale('log')
plt.yscale('log')

plt.title('Time needed per task with IPOT heuristic no preprocessing vs default')
plt.xlabel('Default Time (s)')
plt.ylabel('Time without preprocessing (s)')

plt.xlim(10**-2, 300)
plt.ylim(10**-2, 300)

plt.plot([10**-4, 10**3], [10**-4, 10**3], 'r--', linewidth=1)


import numpy as np
xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()

# Calculate the diagonal lines
log_range = np.logspace(np.log10(xmin), np.log10(xmax), num=50)
factors = np.logspace(np.log10(xmin), np.log10(xmax), num=10)

# Plot diagonal lines above the middle line
for factor in factors:
    plt.plot(log_range, factor * log_range / log_range[0], ls="--", linewidth=0.5, color='gray', alpha=0.5)

# Plot diagonal lines below the middle line
for factor in factors:
    plt.plot(log_range, log_range[0] * log_range / factor, ls="--", linewidth=0.5, color='gray', alpha=0.5)



plt.savefig('scatter_IPOT_no_preprocessing_vs_default.png', dpi=300)

plt.show()


