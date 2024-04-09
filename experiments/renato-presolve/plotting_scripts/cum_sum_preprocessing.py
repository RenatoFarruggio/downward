import json

#print("Hello, world!")

file_path = '2024-02-06_4_no_preprocessing-eval/properties'

algo_name = input("Enter algo name (SEH/DEL/OCP/PHO/IPOT/DPOT): ")

ignored_algorithms = ["config_DEL_default",
                      "config_DEL_no_preprocessing",
                      "config_DPOT_default",
                      "config_DPOT_no_preprocessing",
                      "config_IPOT_default",
                      "config_IPOT_no_preprocessing",
                      "config_OCP_default",
                      "config_OCP_no_preprocessing",
                      "config_PHO_default",
                      "config_PHO_no_preprocessing",
                      "config_SEH_default",
                      "config_SEH_no_preprocessing",
                      ]

default_string = f"config_{algo_name}_default"
no_preprocessing_string = f"config_{algo_name}_no_preprocessing"

if default_string in ignored_algorithms and no_preprocessing_string in ignored_algorithms:
    ignored_algorithms.remove(default_string)
    ignored_algorithms.remove(no_preprocessing_string)
else:
    print(f"ERROR: String {default_string} or {no_preprocessing_string} not found!")
    exit()




success_counter_default = 0
success_counter_no_preprocessing = 0

total_times_default = []
total_times_no_preprocessing = []

with open(file_path, 'r') as file:
    data = json.load(file)



for key in data:
    #print(key)
    error = data[key].get("error", None)
    #print(f"{error=}")
    
    if error == "success":
        algorithm = data[key].get("algorithm", None)
        algorithm = algorithm.split(':')[-1]
        
        if algorithm == default_string:
            success_counter_default += 1
            total_time = data[key].get("total_time", None)
            total_times_default.append(total_time)
        
        elif algorithm == no_preprocessing_string:
            success_counter_no_preprocessing += 1
            total_time = data[key].get("total_time", None)
            total_times_no_preprocessing.append(total_time)
        
        elif algorithm not in ignored_algorithms:
            print("Algorithm not yet handled:", algorithm)
        
        else:
            pass


#        match algorithm:
#            case "config_SEH_default" | "config_SEH_no_preprocessing" \
#                | "config_PHO_default" | "config_PHO_no_preprocessing" \
#                | "config_OCP_default" | "config_OCP_no_preprocessing" \
#                | "config_DEL_default" | "config_DEL_no_preprocessing" \
#                | "config_DPOT_default" | "config_DPOT_no_preprocessing":
#                pass
#            
#            case "config_IPOT_default":
#                success_counter_IPOT_default += 1
#                total_time = data[key].get("total_time", None)
#                total_times_IPOT_default.append(total_time)
#
#            case "config_IPOT_no_preprocessing":
#                success_counter_IPOT_no_preprocessing += 1
#                total_time = data[key].get("total_time", None)
#                total_times_IPOT_no_preprocessing.append(total_time)    
#            
#            case _:
#                print("Algorithm not yet handled:", algorithm)
#                print(data[key].get("total_time", None))

print(f"{success_counter_default=}")
print(f"{success_counter_no_preprocessing=}")


### Plotting ###

import pandas as pd
import matplotlib.pyplot as plt

total_times_default.sort()
total_times_no_preprocessing.sort()


# DataFrame for IPOT default times
df_default = pd.DataFrame({"Time": total_times_default})
df_default['Cumulative Problems Solved default'] = df_default.rank(method='first')

# DataFrame for IPOT no preprocessing times
df_no_preprocessing = pd.DataFrame({"Time": total_times_no_preprocessing})
df_no_preprocessing['Cumulative Problems Solved no preprocessing'] = df_no_preprocessing.rank(method='first')


plt.figure(figsize=(10,6))

plt.plot(df_no_preprocessing['Time'], df_no_preprocessing['Cumulative Problems Solved no preprocessing'], "-x", markevery=32, label=f'{algo_name} no preprocessing')
plt.plot(df_default['Time'], df_default['Cumulative Problems Solved default'], '-D', markevery=32, label=f'{algo_name} default')


plt.title(f'Tasks solved with {algo_name} heuristic with and without preprocessing')
plt.xscale('log')
plt.xlabel('Time (s)')
plt.ylabel('Tasks Solved (Cumulative Sum)')

plt.legend()

plt.savefig(f'cumulative_problems_solved_{algo_name}_no_preprocessing.png', dpi=300)

plt.show()

