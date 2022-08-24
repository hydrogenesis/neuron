import math
import os, sys
import matplotlib.pyplot as plt

def zeros(n):
    ret = range(n)
    return [0 for e in ret]

def dynamicI(t, I):
    if t < 100:
        return I
    elif t > 120 and t < 140:
        return I
    elif t > 160 and t < 180:
        return I
    else:
        return 0

def hh(I,tspan, v, mi, hi, ni,Plot):
    dt = 0.001
    loop = math.ceil(tspan / dt)
    gNa = 120.0
    eNa=115.0
    gK = 36.0
    eK=-12.0
    gL=0.3
    eL=10.6

    t = range(0, loop)
    t = [e * dt for e in t]
    V = zeros(loop)
    m = zeros(loop)
    h = zeros(loop)
    n = zeros(loop)

    V[0]=v
    m[0]=mi
    h[0]=hi
    n[0]=ni
    
    for i in range(loop - 1):
        V[i+1] = V[i] + dt*(gNa*math.pow(m[i],3)*h[i]*(eNa-(V[i]+65)) + gK*math.pow(n[i],4)*(eK-(V[i]+65)) + gL*(eL-(V[i]+65)) + dynamicI(i*dt, I));
        m[i+1] = m[i] + dt*(alphaM(V[i])*(1-m[i]) - betaM(V[i])*m[i]);
        h[i+1] = h[i] + dt*(alphaH(V[i])*(1-h[i]) - betaH(V[i])*h[i]);
        n[i+1] = n[i] + dt*(alphaN(V[i])*(1-n[i]) - betaN(V[i])*n[i]);

    if Plot:
        plt.plot(t, V)
        plt.xlabel('t(ms)')
        plt.ylabel('V(mV)')
        plt.savefig('hh.png')

def alphaM(V):
    return (2.5-0.1*(V+65)) / (math.exp(2.5-0.1*(V+65)) -1)

def betaM(V):
    return 4.0*math.exp(-(V+65)/18)

def alphaH(V):
    return 0.07*math.exp(-(V+65)/20)

def betaH(V):
    return 1.0/(math.exp(3.0-0.1*(V+65))+1)

def alphaN(V):
    return (0.1-0.01*(V+65)) / (math.exp(1-0.1*(V+65)) -1)

def betaN(V):
    return 0.125*math.exp(-(V+65)/80)

if __name__ == '__main__':
    hh(10, 200, -65, 0.5, 0.06, 0.5, True)
