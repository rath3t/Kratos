import matplotlib.pyplot as plt
import re

def read_file(file_name):
    column1 = []
    column2 = []
    point_coordinates = None
    with open(file_name, 'r') as file:
        for line in file:
            if line.startswith('#'):
                extraction = ExtractArrayFromString(line)
                if extraction:
                    point_coordinates = extraction

            else:
                data = line.split()
                column1.append(float(data[0]))
                column2.append(float(data[1]))
    return point_coordinates, [column1, column2]

def ExtractArrayFromString(string):
    pattern = r"\((.*?)\)"
    matches = re.findall(pattern, string)
    floats = []
    for match in matches:
        float_list = re.findall(r"[-+]?\d+\.\d+", match)
        floats.extend([float(num) for num in float_list])
    return floats

file_path_template = 'gid_output/meshXXX/temps.grf'

for i in range(1,3):
    file_path = file_path_template.replace('XXX', str(i))
    point_coordinates, data_pair = read_file(file_path)
    plt.plot(data_pair[0], data_pair[1], label='mesh' + str(i))


plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.title('Temperature evolution at ' + str(point_coordinates))
plt.legend()
plt.grid(True)
plt.savefig('temperature_time_meshes.pdf')