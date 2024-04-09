import json

#print("Hello, world!")

file_path = '2024-03-22_DEL_IP_relaxations-eval/properties'

categories = [#"default",
              #"default_barrier",
              #"use_time_vars",
              #"use_time_vars_barrier",
              "use_integer_vars",
              "use_integer_vars_barrier",
              "use_integer_time_vars",
              "use_integer_time_vars_barrier"]

n = len(categories)

success_counters = [0] * n

total_times = [[] for _ in range(n)]


with open(file_path, 'r') as file:
    data = json.load(file)


for key in data:
    #print(key)
    error = data[key].get("error", None)
    #print(f"{error=}")
    
    if error == "success":
        algorithm = data[key].get("algorithm", None)
        algorithm = algorithm.split(':')[-1]
        

        try:
            i = categories.index(algorithm)

            success_counters[i] += 1
            total_time = data[key].get("total_time", None)
            total_times[i].append(total_time)

        except ValueError:
            print("Algorithm not yet handled:", algorithm)
            #print(data[key].get("total_time", None))
            continue

for i in range(n):
    print(f"Number of {categories[i]}: {success_counters[i]}")
#print(f"{success_counter_default=}")
#print(f"{success_counter_default_barrier=}")
#print(f"{success_counter_use_time_vars=}")


### Plotting ###
import pandas as pd
import matplotlib.pyplot as plt

for item in total_times:
    item.sort()

data_frames = []
for i in range(n):
    data_frame = pd.DataFrame({"Time": total_times[i]})
    data_frame[f'Cumulative Problems Solved {categories[i]}'] = data_frame.rank(method='first')
    data_frames.append(data_frame)

## DataFrame for default times
#df_default = pd.DataFrame({"Time": total_times_default})
#df_default['Cumulative Problems Solved default'] = df_default.rank(method='first')
#
## DataFrame for default times with barrier
#df_default_barrier = pd.DataFrame({"Time": total_times_default_barrier})
#df_default_barrier['Cumulative Problems Solved default barrier'] = df_default_barrier.rank(method='first')
#
## DataFrame for use_time_vars
#df_default_use_time_vars = pd.DataFrame({"Time": total_times_use_time_vars})
#df_default_use_time_vars['Cumulative Problems Solved use_time_vars'] = df_default_use_time_vars.rank(method='first')


plt.figure(figsize=(10,6))

# -^, -D, -o, -*

symbols = ['-^', '--^',
           '-D', '--D',
           '-o', '--o',
           '-*', '--*'
           ]

for i in range(n):
    data_frame = data_frames[i]
    plt.plot(data_frame['Time'], data_frame[f'Cumulative Problems Solved {categories[i]}'], symbols[i], markevery=32, label=categories[i])

#plt.plot(df_default['Time'], df_default['Cumulative Problems Solved default'], '-^', markevery=32, label='default')
#plt.plot(df_default_barrier['Time'], df_default_barrier['Cumulative Problems Solved default barrier'], '-D', markevery=32, label='default barrier')
#plt.plot(df_default_use_time_vars['Time'], df_default_use_time_vars['Cumulative Problems Solved use_time_vars'], '-o', markevery=32, label='use_time_vars')


plt.xscale('log')
plt.xlabel('Time (s)')
plt.ylabel('Problems Solved (Cumulative Sum)')

#plt.grid(True, which="both", ls="--")

plt.legend()

plt.savefig('cumulative_problems_solved_DEL_IP_relaxation.png', dpi=300)

plt.show()

