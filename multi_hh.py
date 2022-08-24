import threading
import math
import time
import os, sys
import matplotlib.pyplot as plt

global cv
global lock

class Params:
    def __init__(self):
        self.cv = threading.Condition()
        self.lock = threading.Lock()
        self.signal = 0
        self.p = 0
        self.total = 0

class Membrane(threading.Thread):
    def __init__(self, params, I, tspan, dt, v, mi, hi, ni, gNa, eNa, gK, eK, gL, eL, g):
        super(Membrane, self).__init__(name="hh")
        self.params = params
        self.i = 0
        self.I = I
        self.dt = dt
        self.loop = math.ceil(tspan / dt)
        self.v = v
        self.mi = mi
        self.hi = hi
        self.ni = ni
        self.gNa = gNa
        self.eNa = eNa
        self.gK = gK
        self.eK = eK
        self.gL = gL
        self.eL = eL
        self.g = g
        self.t = range(0, self.loop)
        self.t = [e * dt for e in self.t]
        self.V = self.zeros(self.loop)
        self.m = self.zeros(self.loop)
        self.h = self.zeros(self.loop)
        self.n = self.zeros(self.loop)
        self.V[0]=self.v
        self.m[0]=self.mi
        self.h[0]=self.hi
        self.n[0]=self.ni
        self.prev = []
        self.next = []

    def dynamicI(self, t, I):
        return I
        if t < 100:
            return I
        elif t > 120 and t < 140:
            return I
        elif t > 160 and t < 180:
            return I
        else:
            return 0

    def zeros(self, n):
        ret = range(n)
        return [0 for e in ret]

    def Append(self, m):
        self.next.append(m)
        m.prev.append(self)

    def run(self):
        while True:
            with self.params.cv:
                self.params.p += 1
                self.params.cv.wait()
            s = self.params.signal
            if s == 0:
                self.Advance()
            elif s == 1:
                break
            else:
                continue

    def Advance(self):
        self.i += 1

    def Advance_Euler(self):
        self.V[self.i+1] = self.V[self.i] + self.dt*(self.gNa*math.pow(self.m[self.i],3)*self.h[self.i]*(self.eNa-(self.V[self.i]+65)) + self.gK*math.pow(self.n[self.i],4)*(self.eK-(self.V[self.i]+65)) + self.gL*(self.eL-(self.V[self.i]+65)) + self.dynamicI(self.i*self.dt, self.I))
        self.m[self.i+1] = self.m[self.i] + self.dt*(self.alphaM(self.V[self.i])*(1-self.m[self.i]) - self.betaM(self.V[self.i])*self.m[self.i])
        self.h[self.i+1] = self.h[self.i] + self.dt*(self.alphaH(self.V[self.i])*(1-self.h[self.i]) - self.betaH(self.V[self.i])*self.h[self.i])
        self.n[self.i+1] = self.n[self.i] + self.dt*(self.alphaN(self.V[self.i])*(1-self.n[self.i]) - self.betaN(self.V[self.i])*self.n[self.i])
        for m in self.next:
            self.V[self.i+1] += self.dt * ((m.V[m.i] - self.V[self.i]) * m.g)
        for m in self.prev:
            self.V[self.i+1] -= self.dt * ((self.V[self.i] - m.V[m.i]) * m.g)
        #print(self.t[self.i+1], self.V[self.i+1], self.m[self.i+1], self.h[self.i+1], self.n[self.i+1])
        #os._exit(0)

    def Advance_Crank_Nicolson(self, iter_num):
        #print("init V:", self.V[self.i+1])
        for it in range(iter_num):
            self.V[self.i+1] = self.V[self.i] + 0.5*self.dt*(self.gNa*math.pow(self.m[self.i+1],3)*self.h[self.i+1]*(self.eNa-(self.V[self.i+1]+65)) + self.gK*math.pow(self.n[self.i+1],4)*(self.eK-(self.V[self.i+1]+65)) + self.gL*(self.eL-(self.V[self.i+1]+65)) + self.dynamicI((self.i+1)*self.dt, self.I) + self.gNa*math.pow(self.m[self.i],3)*self.h[self.i]*(self.eNa-(self.V[self.i]+65)) + self.gK*math.pow(self.n[self.i],4)*(self.eK-(self.V[self.i]+65)) + self.gL*(self.eL-(self.V[self.i]+65)) + self.dynamicI(self.i*self.dt, self.I))
            #v(i+1) = 1/2*(f(i+1) + f(i))*dt + v(i)
            #v[i+1] = 1/2*(gNa*m[i+1]^3*h[i+1]*(eNa-(v[i+1]+65)) + gK*n[i+1]^4*(eK-(v[i+1]+65)) + gL*(eL-(v[i+1]+65)) + I + f(i))*dt + v(i)
            #m[i+1] = 1/2*(alphaM(v[i+1]*(1-m[i+1]) - betaM(v[i+1])*m[i+1] + f(i))*dt + m[i]
            #h[i+1] = 1/2*(alphaH(v[i+1]*(1-h[i+1]) - betaH(v[i+1])*h[i+1] + f(i))*dt + h[i]
            #n[i+1] = 1/2*(alphaN(v[i+1]*(1-n[i+1]) - betaN(v[i+1])*n[i+1] + f(i))*dt + n[i]
            self.m[self.i+1] = self.m[self.i] + 0.5*self.dt*(self.alphaM(self.V[self.i+1])*(1-self.m[self.i+1]) - self.betaM(self.V[self.i+1])*self.m[self.i+1] + self.alphaM(self.V[self.i])*(1-self.m[self.i]) - self.betaM(self.V[self.i])*self.m[self.i])
            self.h[self.i+1] = self.h[self.i] + 0.5*self.dt*(self.alphaH(self.V[self.i+1])*(1-self.h[self.i+1]) - self.betaH(self.V[self.i+1])*self.h[self.i+1] + self.alphaH(self.V[self.i])*(1-self.h[self.i]) - self.betaH(self.V[self.i])*self.h[self.i])
            self.n[self.i+1] = self.n[self.i] + 0.5*self.dt*(self.alphaN(self.V[self.i+1])*(1-self.n[self.i+1]) - self.betaN(self.V[self.i+1])*self.n[self.i+1] + self.alphaN(self.V[self.i])*(1-self.n[self.i]) - self.betaN(self.V[self.i])*self.n[self.i])
            for m in self.next:
                self.V[self.i+1] += self.dt * 0.5 * (((m.V[m.i+1] - self.V[self.i+1]) + (m.V[m.i] - self.V[self.i])) * m.g)
            for m in self.prev:
                self.V[self.i+1] -= self.dt * 0.5 * (((self.V[self.i+1] - m.V[m.i+1]) + (self.V[self.i] - m.V[m.i])) * m.g)
            #print("calc V:", self.V[self.i+1])

    def Plot(self, filename):
        start = 65
        end = 75
        start = 0
        end = 200
        plt.figure(figsize=(10, 10))
        plt.plot(self.t[start*1000:end*1000], self.V[start*1000:end*1000])
        plt.xlabel('t(ms)')
        plt.ylabel('V(mV)')
        plt.savefig(filename)

    def alphaM(self, V):
        return (2.5-0.1*(V+65)) / (math.exp(2.5-0.1*(V+65)) -1)
    
    def betaM(self, V):
        return 4.0*math.exp(-(V+65)/18)
    
    def alphaH(self, V):
        return 0.07*math.exp(-(V+65)/20)
    
    def betaH(self, V):
        return 1.0/(math.exp(3.0-0.1*(V+65))+1)
    
    def alphaN(self, V):
        return (0.1-0.01*(V+65)) / (math.exp(1-0.1*(V+65)) -1)
    
    def betaN(self, V):
        return 0.125*math.exp(-(V+65)/80)

def zeros(n):
    ret = range(n)
    return [0 for e in ret]

def dynamicI(t, I):
    return I
    if t < 100:
        return I
    elif t > 120 and t < 140:
        return I
    elif t > 160 and t < 180:
        return I
    else:
        return 0

def hh(I,tspan, v, mi, hi, ni,Plot):
    params = Params()
    dt = 0.001
    loop = math.ceil(tspan / dt)
    gNa = 120.0
    eNa=115.0
    gK = 36.0
    eK=-12.0
    gL=0.3
    eL=10.6
    g=3
    m = Membrane(params, I, tspan, dt, v, mi, hi, ni, gNa, eNa, gK, eK, gL, eL, g)
    lastM = m
    marray = [m]
    msize = 0
    for i in range(msize):
        m1 = Membrane(params, 0, tspan, dt, v, mi, hi, ni, gNa, eNa, gK, eK, gL, eL, g)
        lastM.Append(m1)
        lastM = m1
        marray.append(m1)

    for i in range(loop - 1):
        for m in marray:
            m.Advance_Euler()
        #for m in marray:
        #    m.Advance_Crank_Nicolson(3)
        for m in marray:
            m.Advance()
    # params.signal = 2
    # params.p = 0
    # params.total = len(marray)
    # for j in marray:
    #     j.start()
    # for i in range(loop - 1):
    #     params.signal = 0
    #     params.p = 0
    #     with params.cv:
    #         params.cv.notify_all()
    #     if i % 1000 == 0: print(i, loop - 1)
    #     while params.p != len(marray):
    #         #params.signal = params.signal
    #         time.sleep(0.00001)

    # params.signal = 1
    # with params.cv:
    #     params.cv.notify_all()
    # for j in marray:
    #     j.join()

    if Plot:
        start = 0
        end = 200
        plt.figure(figsize=(10, 10))
        plt.plot(marray[0].t[start*1000:end*1000], marray[0].V[start*1000:end*1000], marray[-1].t[start*1000:end*1000], marray[-1].V[start*1000:end*1000])
        plt.xlabel('t(ms)')
        plt.ylabel('V(mV)')
        plt.savefig('hh_object.png')

if __name__ == '__main__':
    hh(50, 200, -65, 0.5, 0.06, 0.5, True)
