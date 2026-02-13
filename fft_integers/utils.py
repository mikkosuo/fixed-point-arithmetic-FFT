import numpy as np
import scipy
import sys
from pathlib import Path
import matplotlib.pyplot as plt

import ast

def verify_results(TimeAxis, Input, SNR, EchoPresent):
  # Read C outputs
    with open("output/input_real.txt") as f:
        cont = f.read().strip()
    c_input_real = ast.literal_eval(cont)
    with open("output/input_imag.txt") as f:
        cont = f.read().strip()
    c_input_imag = ast.literal_eval(cont)
    with open("output/matched_filter_real.txt") as f:
        cont = f.read().strip()
    c_matched_filter_real = ast.literal_eval(cont)
    with open("output/matched_filter_imag.txt") as f:
        cont = f.read().strip()
    c_matched_filter_imag = ast.literal_eval(cont)
    with open("output/magnitudes.txt") as f:
        cont = f.read().strip()
    c_magnitudes_after_matched_filter = ast.literal_eval(cont)
    with open("output/output_cfar.txt") as f:
        cont = f.read().strip()
    c_cfar = ast.literal_eval(cont) 

    fig, axs = plt.subplots(3, 2, figsize=(14, 10))  # <- Larger figure
    fig.suptitle(f"SNR ={SNR}dB, EchoPresent={EchoPresent}", fontsize=20)  # large figure-wide title
    fig.subplots_adjust(hspace=0.4, wspace=0.3)  # hspace = vertical spacing, wspace = horizontal spacing
   
   # Python side
    axs[0][0].plot(TimeAxis* 1e6, scale_and_round(Input.real))
    axs[0][0].set_xlabel("Time [µs]")
    axs[0][0].set_title("Input signal (Python)")
    integers_magnitude = np.sqrt(scale_and_round(Input.real)**2 + scale_and_round(Input.imag)**2)
    axs[0][1].plot(TimeAxis* 1e6, integers_magnitude)
    axs[0][1].set_xlabel("Time [µs]")
    axs[0][1].set_title("Input magnitudes (Python)")

    # C
    axs[1][0].plot(TimeAxis* 1e6, c_input_real)
    axs[1][0].set_title("Input signal (C)")
    magnitude = np.sqrt(np.array(c_input_real)**2 + np.array(c_input_imag)**2)
    axs[1][1].plot(TimeAxis* 1e6, magnitude)
    axs[1][1].set_xlabel("Time [µs]")
    axs[1][1].set_title("Input magnitudes (C)")

    axs[2][0].plot(TimeAxis* 1e6, c_matched_filter_real)
    axs[2][0].set_xlabel("Time [µs]")
    axs[2][0].set_title("After matched filter (C)")
    axs[2][1].plot(TimeAxis* 1e6, c_magnitudes_after_matched_filter)
    axs[2][1].plot(TimeAxis* 1e6, c_cfar)
    axs[2][1].set_xlabel("Time [µs]")
    axs[2][1].set_title("Magnitudes and CFAR after matched filter(C)")

    plt.show()


def scale_and_round(x):
    return np.round(x * 2**6).astype(np.int32)

def generate_chirp(P, tau, t, f0, B):
    # Generate unit-amplitude complex chirp
    chirp = scipy.signal.chirp(t,  f0=f0, t1=tau, f1=f0 + B, method='linear',
        phi=-90,
        complex=True
    )
    # Normalize to unit average power
    chirp /= np.sqrt(np.mean(np.abs(chirp)**2))
    # Scale to desired power
    chirp *= np.sqrt(P)
    return chirp

def generate_echo_signal(P, tau, t, f0, B, t_delay, SNR_dB, fs, N_total, N_chirp, EchoPresent = True):
    chirp_back = generate_chirp(P, tau, t, f0, B)
    rx =  np.zeros(N_total, dtype=np.complex64)
    noise_power = P / (10**(SNR_dB / 10))
    noise_std = np.sqrt(noise_power / 2)
    noise = noise_std * (np.random.randn(N_total) + 1j * np.random.randn(N_total))
    if EchoPresent:
        delayIndex = int(fs * t_delay) 
        rx[delayIndex:N_chirp+delayIndex] = chirp_back
    rx = rx + noise
    return np.array(rx)

def cfar_thresholds(Input, GAP, REF, bias):
  N = Input.size
  output = np.zeros(N)
  for center_index in range(GAP + REF, N - (GAP + REF)):
    min_index = center_index - (GAP + REF)
    min_guard = center_index - GAP
    max_index = center_index + GAP + REF + 1
    max_guard = center_index + GAP + 1

    lower = Input[min_index:min_guard]
    upper = Input[max_guard:max_index]

    mean = np.mean(np.concatenate((lower,upper)))
    output[center_index] = mean * bias
  return output

def matched_filter(ReceivedSignal, Transmitted):
  h = np.conj(Transmitted[::-1])
  return np.convolve(ReceivedSignal, h, mode="same")