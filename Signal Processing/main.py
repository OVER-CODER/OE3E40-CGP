import numpy as np
import matplotlib.pyplot as plt

# --- Parameters ---
fs = 1000       # Sampling frequency (Hz)
f1 = 50         # Frequency 1 (Hz)
f2 = 120        # Frequency 2 (Hz)
N1 = 256        # FFT length 1
N2 = 1024       # FFT length 2

# --- Generation of input sequence ---
t = np.arange(0, 1, 1/fs)
x = np.sin(2*np.pi*f1*t) + 2*np.sin(2*np.pi*f2*t) - np.random.randn(len(t))

# --- Generation of PSD for two different FFT lengths ---
Pxx1 = np.abs(np.fft.fft(x, N1))**2 / (N1 + 1)
Pxx2 = np.abs(np.fft.fft(x, N2))**2 / (N2 + 1)

# --- Frequency vectors ---
f1_axis = np.arange(0, N1) / N1 * fs
f2_axis = np.arange(0, N2) / N2 * fs

# --- Plot the PSD ---
plt.figure(figsize=(10, 6))

plt.subplot(2, 1, 1)
plt.plot(f1_axis, 10*np.log10(Pxx1))
plt.grid(True)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power Spectrum (dB)')
plt.title(f'PSD with FFT length {N1}')

plt.subplot(2, 1, 2)
plt.plot(f2_axis, 10*np.log10(Pxx2))
plt.grid(True)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power Spectrum (dB)')
plt.title(f'PSD with FFT length {N2}')

plt.tight_layout()
plt.show()
