from utils import *

Q = 15
SCALE = 32767
N = 1024 
fs = 1024

def multi_frequence_example():
  t = np.arange(N)
  # --- Frequencies ---
  f1 = 300          # strong tone
  f2 = 123          # medium tone
  f3 = 73.5         # non-integer bin (leakage test)
  f4 = 50           # AM carrier
  f_mod = 7         # AM modulation freq
  
  # --- Build components ---
  
  # Strong complex tone
  sig1 = 0.7 * np.exp(1j * 2*np.pi*f1*t/fs)
  
  # Medium complex tone with phase offset
  sig2 = 0.4 * np.exp(1j * (2*np.pi*f2*t/fs + np.pi/4))
  
  # Non-bin aligned tone (causes leakage)
  sig3 = 0.25 * np.exp(1j * 2*np.pi*f3*t/fs)
  
  # Amplitude-modulated complex tone
  carrier = np.exp(1j * 2*np.pi*f4*t/fs)
  mod = (1 + 0.5*np.sin(2*np.pi*f_mod*t/fs))
  sig4 = 0.3 * mod * carrier
  
  # Add small noise
  np.random.seed(0)
  noise = 0.05 * (np.random.randn(N) + 1j*np.random.randn(N))
  
  # Combine everything
  x_complex = sig1 + sig2 + sig3 + sig4 + noise
  
  # Normalize to avoid clipping
  x_complex /= np.max(np.abs(x_complex))
  x = np.zeros((N,2), dtype=np.int32)
  x[:,0] = np.round(np.real(x_complex) * SCALE).astype(np.int32)
  x[:,1] = np.round(np.imag(x_complex) * SCALE).astype(np.int32)
  return x

def single_frequence_example():
  f = 300 
  t = np.arange(N) 
  x_real = np.cos(2*np.pi*f*t/fs)
  x_imag = np.sin(2*np.pi*f*t/fs)
  x = np.zeros((N,2), dtype=np.int32)
  x[:,0] = (x_real * SCALE).astype(np.int32)
  x[:,1] = (x_imag * SCALE).astype(np.int32)
  return x


def main():
  x = multi_frequence_example()
  # Write the header
  base_dir = Path(__file__).resolve().parent
  output_file = base_dir / 'test_data.h'
  
  with open(output_file, 'w') as fd:
    fd.write("#include <inttypes.h>\n");

    fd.write(f"int16_t signal[{N}][2] = {{\n")
    for i in range(N):
      fd.write(f"{{ {x[i][0]:g}, {x[i][1]:g} }},\n")
    fd.write("};\n")
    fd.write("\n")
  
  print("Done")
if __name__ == "__main__":
  main()            


