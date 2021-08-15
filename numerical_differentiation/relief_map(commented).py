import struct  # for reading file
import math  # for calculating intensity
import matplotlib.pyplot as plt  # for plotting

phi = math.pi / 6  # angle of the incident light
h = 83  # distance between grid points in meters
grid_dim = 1201  # dimensions of the file data
file_name = 'N46E006.hgt'  # file name of the map data as a str

f = open(file_name, 'rb')

# create 2D array of the elevation at each data point
elevation_map = []  # initialize row array
for y in range(grid_dim):  # iterate through each row of data
    elevation_row = []  # initialize column array
    for x in range(grid_dim):  # iterate through each column of data
        elevation_row.append(struct.unpack('>h', f.read(2))[0])  # add elevation data to column array
    elevation_map.append(elevation_row)  # add column array to row array

# create a 2D array of the gradient vector and a 2D array of the illumination intensity
gradient_map = []  # initialize row array
intensity_map = []  # initialize row array
for y in range(grid_dim):  # iterate through each row of data
    gradient_row = []  # initialize column array
    intensity_row = []  # initialize column array
    for x in range(grid_dim):  # iterate through each column of data
        gradient = []  # initialize gradient vector
        if x != grid_dim - 1:  # if the data point is not an edge case on the right side of the map
            dw_dx = (elevation_map[y][x + 1] - elevation_map[y][x]) / h  # calculate the rate of altitude change towards the positive x direction
        else:  # if the data point is an edge case on the right side of the map
            dw_dx = (elevation_map[y][x] - elevation_map[y][x - 1]) / h  # calculate the rate of altitude change from the negative x direction
        if y != grid_dim - 1:  # if the data point is not an edge case on the bottom of the map
            dw_dy = (elevation_map[y][x] - elevation_map[y + 1][x]) / h  # calculate the rate of altitude change from the negative y direction
        else:  # if the data point is an edge case on the bottom of the map
            dw_dy = (elevation_map[y - 1][x] - elevation_map[y][x]) / h  # calculate the rate of altitude change towards the positive y direction
        gradient.append(dw_dx)  # add dw/dx as the first term of the gradient vector
        gradient.append(dw_dy)  # add dw/dy as the second term of the gradient vector
        gradient.append(-1)  # add -1 as the third term of the gradient vector
        gradient_row.append(gradient)  # add gradient vector to the column array

        # calculate the illumination intensity from the gradient vector
        dividend = math.cos(phi) * (gradient[0]) + math.sin(phi) * (gradient[1])
        divisor = math.sqrt(gradient[0] ** 2 + gradient[1] ** 2 + 1)
        intensity_row.append(dividend / divisor)  # add intensity value to the column array
    gradient_map.append(gradient_row)  # add column array to the row array
    intensity_map.append(intensity_row)  # add column array to the row array

x_ticks = []  # initialize array for the original axis ticks
x_ticks_new = []  # initialize array for the new axis ticks
y_ticks = []  # initialize array for the original axis ticks
y_ticks_new = []  # initialize array for the new axis ticks

for i in range(0, grid_dim, int((grid_dim - 1) / 10)):  # iterate through every tenth of the map data's dimension
    x_ticks.append(i)  # add x coordinate to original axis array
    y_ticks.append(grid_dim - i)  # add y coordinate, adjusted to positive y-axis, to original axis array
    x_ticks_new.append(int(file_name[4:7]) + i / (grid_dim - 1))  # calculate longitude at the corresponding x coordinate using the file name and add longitude to new axis array
    y_ticks_new.append(int(file_name[1:3]) + i / (grid_dim - 1))  # calculate latitude at the corresponding y coordinate using the file name and add latitude to new axis array

plt.figure()
plt.imshow(elevation_map, 'gray', vmin=0)  # plot elevation map and use vmin to dismiss broken data
plt.title("Elevation map")
plt.xticks(x_ticks, x_ticks_new)
plt.yticks(y_ticks, y_ticks_new)
plt.xlabel("Latitude E($^\circ$)")
plt.ylabel("Longitude N($^\circ$)")
plt.colorbar().set_label('Elevation (m)')
plt.savefig("Elevation map")

plt.figure()
plt.imshow(intensity_map, 'gray', vmin=-0.5, vmax=0.5)  # plot elevation map and use vmin / vmax to dismiss broken data
plt.title("Intensity map of surface illumination")
plt.xticks(x_ticks, x_ticks_new)
plt.yticks(y_ticks, y_ticks_new)
plt.xlabel("Latitude E($^\circ$)")
plt.ylabel("Longitude N($^\circ$)")
plt.colorbar().set_label('Intensity')
plt.savefig("Intensity map of surface illumination")

plt.show()
