elevation_map = [
    [1, 1, 1, 1, 1],
    [1, 0, 0, 0, 1],
    [1, 0, 0, 0, 1],
    [1, 0, 0, 0, 1],
    [1, 1, 1, 1, 1]
]
grid_dim = 5  # dimensions of the file data
h = 1

dw_dx_map = []
dw_dy_map = []
for y in range(grid_dim):
    dw_dx_row = []
    dw_dy_row = []
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
        dw_dx_row.append(dw_dx)
        dw_dy_row.append(dw_dy)
    dw_dx_map.append(dw_dx_row)
    dw_dy_map.append(dw_dy_row)

print(dw_dx_map)
print(dw_dy_map)
