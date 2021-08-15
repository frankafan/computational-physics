import matplotlib.pyplot as plt
import math
import numpy as np

FILE_NAME = 'blur.txt'
sigma = 25
error = 10e-3

image_matrix = []
with open(FILE_NAME, 'r') as image_file:
    line = image_file.readline()
    while line != '':
        line_lst = list(map(float, line.strip().split(' ')))
        image_matrix.append(line_lst)
        line = image_file.readline()
image_height = len(image_matrix)
image_width = len(image_matrix[0])


def f(x, y):
    return math.exp(-(x ** 2 + y ** 2) / (2 * sigma ** 2))


gaussian_matrix = []
for y in range(image_height):
    yp = y
    if yp > image_height / 2:
        yp -= image_height
    gaussian_row = []
    for x in range(image_width):
        xp = x
        if xp > image_width / 2:
            xp -= image_width
        gaussian_row.append(f(xp, yp))
    gaussian_matrix.append(gaussian_row)

fourier_deconvoluted = []
fourier_image = np.fft.rfft2(image_matrix)
fourier_gaussian = np.fft.rfft2(gaussian_matrix)

for y in range(len(fourier_image)):
    deconvoluted_row = []
    for x in range(len(fourier_image[0])):
        if fourier_gaussian[y][x] >= error:
            deconvoluted_row.append(fourier_image[y][x] / fourier_gaussian[y][x])
        else:
            deconvoluted_row.append(fourier_image[y][x])
    fourier_deconvoluted.append(deconvoluted_row)

deconvoluted_matrix = np.fft.irfft2(fourier_deconvoluted)

plt.figure()
plt.imshow(image_matrix, 'gray')
plt.title("Image file shown with density plot")
plt.colorbar().set_label('Brightness')
plt.savefig("Image file shown with density plot")

plt.figure()
plt.imshow(gaussian_matrix, 'gray')
plt.title("Point spread function visualization")
plt.colorbar().set_label('Gaussian')
plt.savefig("Point spread function visualization")

plt.figure()
plt.imshow(deconvoluted_matrix, 'gray')
plt.title("Deconvoluted image")
plt.colorbar().set_label('Brightness')
plt.savefig("Deconvoluted image")

plt.show()
