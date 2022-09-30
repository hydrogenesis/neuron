#!/usr/bin/python3
import struct
import os, sys

def hex_to_double(f):
    return struct.unpack('!d', bytes.fromhex(f.lstrip("0x")))[0]

print(hex_to_double(sys.argv[1]))
