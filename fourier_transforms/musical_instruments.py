import numpy as np
import matplotlib.pyplot as plt

FILE_NAME = 'piano.txt'
FREQUENCY = 44100

waveform_lst = []
with open(FILE_NAME, 'r') as waveform_file:
    line = waveform_file.readline()
    while line != '':
        waveform_lst.append(float(line.strip()))
        line = waveform_file.readline()
time_lst = []
for i in range(len(waveform_lst)):
    time_lst.append(i / FREQUENCY)

fourier_coefficients = np.fft.rfft(waveform_lst)
fourier_frequency = np.fft.rfftfreq(len(waveform_lst))
frequency_lst = []
for i in range(len(fourier_frequency)):
    frequency_lst.append(abs(fourier_frequency[i] * FREQUENCY))

plt.figure()
plt.plot(time_lst, waveform_lst)
plt.xlabel('Time (s)')
plt.title("Waveform plot of " + FILE_NAME[:-4])
plt.savefig("Waveform plot of " + FILE_NAME[:-4])

plt.figure()
plt.plot(fourier_coefficients)
plt.xlabel('k')
plt.ylabel('|$c_k$|')
plt.title("Fourier transformed waveform plot of " + FILE_NAME[:-4])
plt.savefig("Fourier transformed waveform plot of " + FILE_NAME[:-4])

plt.figure()
plt.plot(frequency_lst, fourier_coefficients)
plt.xlabel('Frequency (Hz)')
plt.ylabel('|$c_k$|')
plt.title("Fourier coefficients with respect to frequency " + FILE_NAME[:-4])
plt.savefig("Fourier coefficients with respect to frequency " + FILE_NAME[:-4])

plt.show()
