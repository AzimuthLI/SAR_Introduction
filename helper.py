# Import necessary package here.
import numpy as np
import struct



def open_dat_file(path: str) -> np.ndarray:
    print(f"Opening file: {path}")
    with open(path, 'rb') as rfile:
        rng = struct.unpack('>l', rfile.read(4))[0]
        azm = struct.unpack('>l', rfile.read(4))[0]
        print(f" - Range:   {rng}")
        print(f" - Azimuth: {azm}")
        dt = np.dtype('>c8')
        arr = np.fromfile(rfile, dtype=dt)
        arr = arr.astype(np.complex_)
    return arr.reshape((azm, rng))