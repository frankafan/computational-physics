import struct
import math
import matplotlib.pyplot as plt

phi = math.pi / 6  # angle of the incident light
h = 83  # distance between grid points in meters
grid_dim = 1201  # dimensions of the file data
file_name = 'N46E006.hgt'  # file name of the map data as a str

f = open(file_name, 'rb')

elevation_map = []
for y in range(grid_dim):
    elevation_row = []
    for x in range(grid_dim):
        elevation_row.append(struct.unpack('>h', f.read(2))[0])
    elevation_map.append(elevation_row)

gradient_map = []
intensity_map = []
for y in range(grid_dim):
    gradient_row = []
    intensity_row = []
    for x in range(grid_dim):
        gradient = []
        if x != grid_dim - 1:
            dw_dx = (elevation_map[y][x + 1] - elevation_map[y][x]) / h
        else:
            dw_dx = (elevation_map[y][x] - elevation_map[y][x - 1]) / h
        if y != grid_dim - 1:
            dw_dy = (elevation_map[y][x] - elevation_map[y + 1][x]) / h
        else:
            dw_dy = (elevation_map[y - 1][x] - elevation_map[y][x]) / h
        gradient.append(dw_dx)
        gradient.append(dw_dy)
        gradient.append(-1)
        gradient_row.append(gradient)

        dividend = math.cos(phi) * (gradient[0]) + math.sin(phi) * (gradient[1])
        divisor = math.sqrt(gradient[0] ** 2 + gradient[1] ** 2 + 1)
        intensity_row.append(dividend / divisor)
    gradient_map.append(gradient_row)
    intensity_map.append(intensity_row)

x_ticks = []
x_ticks_new = []
y_ticks = []
y_ticks_new = []

for i in range(0, grid_dim, int((grid_dim - 1) / 10)):
    x_ticks.append(i)
    y_ticks.append(grid_dim - i)
    x_ticks_new.append(int(file_name[4:7]) + i / (grid_dim - 1))
    y_ticks_new.append(int(file_name[1:3]) + i / (grid_dim - 1))

plt.figure()
plt.imshow(elevation_map, 'gray', vmin=0)
plt.title("Elevation map")
plt.xticks(x_ticks, x_ticks_new)
plt.yticks(y_ticks, y_ticks_new)
plt.xlabel("Latitude E($^\circ$)")
plt.ylabel("Longitude N($^\circ$)")
plt.colorbar().set_label('Elevation (m)')
plt.savefig("Elevation map")

plt.figure()
plt.imshow(intensity_map, 'gray', vmin=-0.5, vmax=0.5)
plt.title("Intensity map of surface illumination")
plt.xticks(x_ticks, x_ticks_new)
plt.yticks(y_ticks, y_ticks_new)
plt.xlabel("Latitude E($^\circ$)")
plt.ylabel("Longitude N($^\circ$)")
plt.colorbar().set_label('Intensity')
plt.savefig("Intensity map of surface illumination")

plt.show()
