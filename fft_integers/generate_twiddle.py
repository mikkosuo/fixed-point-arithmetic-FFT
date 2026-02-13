import math
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

Q = 15
SCALE = 32767
def generate_twiddle_fft(N):
    lut = np.zeros((N//2, 2), dtype=np.int32)
    for k in range(N//2):
        angle = -2 * math.pi * k / N
        lut[k,0] = int(math.cos(angle) * SCALE)
        lut[k,1] = int(math.sin(angle) * SCALE)
    return lut

def main(N):
    
    N = int(N)
    twiddle = generate_twiddle_fft(N)

    twiddle_shape_0, twiddle_shape_1 = twiddle.shape
    # Write the header
    base_dir = Path(__file__).resolve().parent
    output_file = base_dir / "twiddle.h"

    with open(output_file, "w") as fd:
        fd.write("#include <stdint.h>\n\n")

        fd.write(f"#define FFT_N {N}\n\n")

        fd.write(f"static const int16_t twiddle[{N//2}][2] = {{\n")
        for k in range(N//2):
            fd.write(f"{{ {twiddle[k, 0]}, {twiddle[k, 1]} }},\n")
        fd.write("};\n\n")

    print(f"Generated twiddle.h for N = {N}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 gen_twiddle.py <N>")
        sys.exit(1)
    main(sys.argv[1])    


