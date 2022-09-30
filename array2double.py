#!/usr/bin/python3
import struct
import os, sys
import matplotlib.pyplot as plt

def hex_to_double(f):
    return struct.unpack('!d', bytes.fromhex(f.lstrip("0x")))[0]

bytes_array = sys.argv[1]
values = bytes_array.split(',')
result = []
for value in values:
    result.append(hex_to_double(value))
plt.plot(result)
plt.ylabel('V(mV)')
plt.savefig('sol.png')
