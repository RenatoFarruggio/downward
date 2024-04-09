import json

#print("Hello, world!")

file_path = '2024-03-13_OCP_barrier-eval/properties'

with open(file_path, 'r') as file:
    data = json.load(file)

success_counter_default = 0
success_counter_barrier = 0

total_times_default = []
total_times_barrier = []


for key in data:
    #print(key)
    error = data[key].get("error", None)
    #print(f"{error=}")
    
    if error == "success":
        algorithm = data[key].get("algorithm", None)
        algorithm = algorithm.split(':')[-1]
        # Only consider SEH for now
        if algorithm == "barrier":
            #print(f"{algorithm=}")
            success_counter_barrier += 1
            total_time = data[key].get("total_time", None)
            #print(f"{total_time=}")
            total_times_barrier.append(total_time)
        
        elif algorithm == "default":
            success_counter_default += 1
            total_time = data[key].get("total_time", None)
            total_times_default.append(total_time)
        
        else:
            print("Algorithm not yet handled:", algorithm)

    #print()

print(f"{success_counter_barrier=}")
print(f"{success_counter_default=}")


### Plotting ###
# I have a list of total_times by now (only using SEH for now)

import pandas as pd
import matplotlib.pyplot as plt

total_times_barrier.sort()
total_times_default.sort()


# DataFrame for OCP barrier times
df_barrier = pd.DataFrame({"Time": total_times_barrier})
df_barrier['Cumulative Problems Solved OCP barrier'] = df_barrier.rank(method='first')

# DataFrame for OCP default times
df_default = pd.DataFrame({"Time": total_times_default})
df_default['Cumulative Problems Solved OCP default'] = df_default.rank(method='first')


plt.figure(figsize=(10,6))

plt.plot(df_barrier['Time'], df_barrier['Cumulative Problems Solved OCP barrier'], "-x", markevery=32, label='barrier')
plt.plot(df_default['Time'], df_default['Cumulative Problems Solved OCP default'], "-D", markevery=32, label='default')


plt.title('Tasks solved with OCP heuristic with default and barrier algorithm')
plt.xscale('log')
plt.xlabel('Time (s)')
plt.ylabel('Tasks Solved (Cumulative Sum)')

#plt.grid(True, which="both", ls="--")

plt.legend()

plt.savefig('cumulative_problems_solved_OCP_barrier.png', dpi=300)

plt.show()

