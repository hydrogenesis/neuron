#!/usr/bin/python3
import struct
import os, sys

def double_to_hex(f):
    return hex(struct.unpack('<Q', struct.pack('<d', f))[0]).rstrip('L')

print(double_to_hex(float(sys.argv[1])))
