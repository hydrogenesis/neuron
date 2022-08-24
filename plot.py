import os, sys 
import matplotlib.pyplot as plt 

if len(sys.argv) < 2:
    filename = "hh.txt"
else:
    filename = sys.argv[1]
f = open(filename, "r")
T1 = []
V1 = []
T2 = []
V2 = []
for line in f:
    t1, v1, t2, v2 = line.strip().split("\t")
    T1.append(float(t1))
    V1.append(float(v1))
    T2.append(float(t2))
    V2.append(float(v2))

plt.figure(figsize = (10, 10))
plt.plot(T1, V1, T2, V2)
plt.xlabel('t(ms)')
plt.ylabel('V(mV)')
plt.savefig('hh_cpp.png')
