import matplotlib.pyplot as plt

a = ['wildtype (not infected)']
b = ['wildtype (infected)']
c = ['ibr5 knockout']
d = ['R gene knockout']
num_bins = 3

n, bins, patches = plt.hist(a, num_bins, facecolor='blue', alpha=1, ec='black')
n, bins, patches = plt.hist(b, num_bins, facecolor='yellow', alpha=1,
                            ec='black')
n, bins, patches = plt.hist(c, num_bins, facecolor='green', alpha=1, ec='black')
n, bins, patches = plt.hist(d, num_bins, facecolor='red', alpha=1, ec='black')

plt.title(
    "Expected observation before pathogen resistance assay")
plt.xlabel("Plant type")
plt.ylabel("Number of pathogen colonies")
plt.ylim(0,4)
plt.show()
