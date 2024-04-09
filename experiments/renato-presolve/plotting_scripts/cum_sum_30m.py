import json

#print("Hello, world!")

file_path = '2024-04-02_3_baseline_methods_used_long-eval/properties'

with open(file_path, 'r') as file:
    data = json.load(file)

success_counter_seh = 0
success_counter_del = 0
success_counter_ocp = 0
success_counter_pho = 0
success_counter_ipot = 0
success_counter_dpot = 0

total_times_seh = []
total_times_del = []
total_times_ocp = []
total_times_pho = []
total_times_ipot = []
total_times_dpot = []


for key in data:
    #print(key)
    error = data[key].get("error", None)
    #print(f"{error=}")
    
    if error == "success":
        algorithm = data[key].get("algorithm", None)
        algorithm = algorithm.split(':')[-1]
        # Only consider SEH for now
        if algorithm == "SEH":
            #print(f"{algorithm=}")
            success_counter_seh += 1
            total_time = data[key].get("total_time", None)
            #print(f"{total_time=}")
            total_times_seh.append(total_time)
        
        elif algorithm == "DEL":
            success_counter_del += 1
            total_time = data[key].get("total_time", None)
            total_times_del.append(total_time)
        
        elif algorithm == "OCP":
            success_counter_ocp += 1
            total_time = data[key].get("total_time", None)
            total_times_ocp.append(total_time)
        
        elif algorithm == "PHO":
            success_counter_pho += 1
            total_time = data[key].get("total_time", None)
            total_times_pho.append(total_time)

        elif algorithm == "IPOT":
            success_counter_ipot += 1
            total_time = data[key].get("total_time", None)
            total_times_ipot.append(total_time)

        elif algorithm == "DPOT":
            success_counter_dpot += 1
            total_time = data[key].get("total_time", None)
            total_times_dpot.append(total_time)

        else:
            print("Algorithm not yet handled:", algorithm)

    #print()

print(f"{success_counter_seh=}")
print(f"{success_counter_del=}")
print(f"{success_counter_ocp=}")
print(f"{success_counter_pho=}")
print(f"{success_counter_ipot=}")
print(f"{success_counter_dpot=}")


### Plotting ###
# I have a list of total_times by now (only using SEH for now)

import pandas as pd
import matplotlib.pyplot as plt

total_times_seh.sort()
total_times_del.sort()
total_times_ocp.sort()
total_times_pho.sort()
total_times_ipot.sort()
total_times_dpot.sort()


# DataFrame for SEH times
df_seh = pd.DataFrame({"Time": total_times_seh})
df_seh['Cumulative Problems Solved SEH'] = df_seh.rank(method='first')

# DataFrame for DEL times
df_del = pd.DataFrame({"Time": total_times_del})
df_del['Cumulative Problems Solved DEL'] = df_del.rank(method='first')

# DataFrame for OCP times
df_ocp = pd.DataFrame({"Time": total_times_ocp})
df_ocp['Cumulative Problems Solved OCP'] = df_ocp.rank(method='first')

# DataFrame for PHO times
df_pho = pd.DataFrame({"Time": total_times_pho})
df_pho['Cumulative Problems Solved PHO'] = df_pho.rank(method='first')

# DataFrame for IPOT times
df_ipot = pd.DataFrame({"Time": total_times_ipot})
df_ipot['Cumulative Problems Solved IPOT'] = df_ipot.rank(method='first')

# DataFrame for DPOT times
df_dpot = pd.DataFrame({"Time": total_times_dpot})
df_dpot['Cumulative Problems Solved DPOT'] = df_dpot.rank(method='first')

plt.figure(figsize=(10,6))

plt.plot(df_ipot['Time'], df_ipot['Cumulative Problems Solved IPOT'], "-D", markevery=32, label='IPOT')
plt.plot(df_dpot['Time'], df_dpot['Cumulative Problems Solved DPOT'], "-^", markevery=32, label='DPOT')
plt.plot(df_ocp['Time'], df_ocp['Cumulative Problems Solved OCP'], "-p", markevery=32, label='OCP')
plt.plot(df_seh['Time'], df_seh['Cumulative Problems Solved SEH'], "-x", markevery=32, label='SEH')
plt.plot(df_pho['Time'], df_pho['Cumulative Problems Solved PHO'], "-*", markevery=32, label='PHO')
plt.plot(df_del['Time'], df_del['Cumulative Problems Solved DEL'], "-o", markevery=32, label='DEL')


plt.title('Tasks solved across all heuristics with default settings and max runtime 30m')
plt.xscale('log')
plt.xlabel('Time (s)')
plt.ylabel('Tasks Solved (Cumulative Sum)')

#plt.grid(True, which="both", ls="--")

plt.legend()

plt.savefig('cumulative_problems_solved_baseline_30m.png', dpi=300)

plt.show()

