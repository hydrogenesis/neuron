import math
import os, sys
import matplotlib.pyplot as plt

def zeros(n):
    ret = range(n)
    return [0 for e in ret]

def dynamicI(t, I):
    if 1 < t and t < 1.1:
        return I
    else:
        return 0
    # if t < 100:
    #     return I
    # elif t > 120 and t < 140:
    #     return I
    # elif t > 160 and t < 180:
    #     return I
    # else:
    #     return 0

def hh(I,tspan, v, mi, hi, ni, celsius, Plot):
    area = 4*3.1415927*10*10
    dt = 0.001
    loop = math.ceil(tspan / dt)
    gNa = 0.12 * area
    eNa=115.0
    gK = 0.036 * area
    eK=-12.0
    gL=0.0003 * area
    eL=10.7
    q10 = math.pow(3, ((celsius - 6.3)/10))

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
        m[i+1] = m[i] + dt*q10*(alphaM(V[i])*(1-m[i]) - betaM(V[i])*m[i]);
        h[i+1] = h[i] + dt*q10*(alphaH(V[i])*(1-h[i]) - betaH(V[i])*h[i]);
        n[i+1] = n[i] + dt*q10*(alphaN(V[i])*(1-n[i]) - betaN(V[i])*n[i]);

    if Plot:
        plt.plot(t, V)
        plt.xlabel('t(ms)')
        plt.ylabel('V(mV)')
        plt.savefig('hh.png')

def vtrap(x, y):
    if abs(x/y) < 1e-6:
        return y*(1 - x/y/2.0)
    else:
        return x/(math.exp(x/y) - 1)

def alphaM(V):
    # return (2.5-0.1*(V+65)) / (math.exp(2.5-0.1*(V+65)) -1)
    return .1 * vtrap(-(V+40),10)

def betaM(V):
    return 4.0*math.exp(-(V+65)/18)

def alphaH(V):
    return 0.07*math.exp(-(V+65)/20)

def betaH(V):
    return 1.0/(math.exp(3.0-0.1*(V+65))+1)

def alphaN(V):
    #return (0.1-0.01*(V+65)) / (math.exp(1-0.1*(V+65)) -1)
    return .01*vtrap(-(V+55),10)

def betaN(V):
    return 0.125*math.exp(-(V+65)/80)

if __name__ == '__main__':
    V = -65
    minf = alphaM(V)/(alphaM(V)+betaM(V))
    hinf = alphaH(V)/(alphaH(V)+betaH(V))
    ninf = alphaN(V)/(alphaN(V)+betaN(V))
    print(minf, hinf, ninf)
    hh(100000, 10, V, minf, hinf, ninf, 6.3, True)
