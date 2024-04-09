import json

#print("Hello, world!")

file_path = '2024-02-06_4_no_preprocessing-eval/properties'


success_counter_IPOT_default = 0
success_counter_IPOT_no_preprocessing = 0

total_times_IPOT_default = []
total_times_IPOT_no_preprocessing = []

with open(file_path, 'r') as file:
    data = json.load(file)



for key in data:
    #print(key)
    error = data[key].get("error", None)
    #print(f"{error=}")
    
    if error == "success":
        algorithm = data[key].get("algorithm", None)
        algorithm = algorithm.split(':')[-1]
        
        match algorithm:
            case "config_SEH_default" | "config_SEH_no_preprocessing" \
                | "config_PHO_default" | "config_PHO_no_preprocessing" \
                | "config_OCP_default" | "config_OCP_no_preprocessing" \
                | "config_DEL_default" | "config_DEL_no_preprocessing" \
                | "config_DPOT_default" | "config_DPOT_no_preprocessing":
                pass
            
            case "config_IPOT_default":
                success_counter_IPOT_default += 1
                total_time = data[key].get("total_time", None)
                total_times_IPOT_default.append(total_time)

            case "config_IPOT_no_preprocessing":
                success_counter_IPOT_no_preprocessing += 1
                total_time = data[key].get("total_time", None)
                total_times_IPOT_no_preprocessing.append(total_time)    
            
            case _:
                print("Algorithm not yet handled:", algorithm)
                print(data[key].get("total_time", None))

print(f"{success_counter_IPOT_default=}")
print(f"{success_counter_IPOT_no_preprocessing=}")


### Plotting ###

import pandas as pd
import matplotlib.pyplot as plt

total_times_IPOT_default.sort()
total_times_IPOT_no_preprocessing.sort()


# DataFrame for IPOT default times
df_IPOT_default = pd.DataFrame({"Time": total_times_IPOT_default})
df_IPOT_default['Cumulative Problems Solved IPOT default'] = df_IPOT_default.rank(method='first')

# DataFrame for IPOT no preprocessing times
df_IPOT_no_preprocessing = pd.DataFrame({"Time": total_times_IPOT_no_preprocessing})
df_IPOT_no_preprocessing['Cumulative Problems Solved IPOT no preprocessing'] = df_IPOT_no_preprocessing.rank(method='first')


plt.figure(figsize=(10,6))

plt.plot(df_IPOT_no_preprocessing['Time'], df_IPOT_no_preprocessing['Cumulative Problems Solved IPOT no preprocessing'], "-x", markevery=32, label='IPOT no preprocessing')
plt.plot(df_IPOT_default['Time'], df_IPOT_default['Cumulative Problems Solved IPOT default'], '-D', markevery=32, label='IPOT default')


plt.title('Tasks solved with IPOT heuristic with and without preprocessing')
plt.xscale('log')
plt.xlabel('Time (s)')
plt.ylabel('Problems Solved (Cumulative Sum)')

plt.legend()

plt.savefig('cumulative_problems_solved_IPOT_no_preprocessing.png', dpi=300)

plt.show()

