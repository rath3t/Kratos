import matplotlib.pyplot as plt
import os

def read_file(file_name):
    column1 = []
    column2 = []
    with open(file_name, 'r') as file:
        for line in file:
            if line.startswith('#'):
                pass
            else:
                data = line.split()
                column1.append(float(data[0]))
                column2.append(float(data[1]))
    return column1, column2

file1 = "gid_output/mesh1/temperature_time.grf"
column1_1, column2_1 = read_file(file1)

file2 = "gid_output/mesh2/temperature_time.grf"
column1_2, column2_2 = read_file(file2)

plt.plot(column1_1, column2_1, label='mesh 1')
plt.plot(column1_2, column2_2, label='mesh 2')

plt.xlabel('Time (s)')
plt.ylabel('Temperature (ºC)')
plt.title('Temperature - time at (-0.009,-0.009,0.005)')
plt.legend()
plt.grid(True)
plt.savefig('temperature_time_meshes.pdf')