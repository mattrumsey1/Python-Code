#This code has been used for all orders of the Runge-Kutta method so parts are 
#hashtagged out depending on what graph is being created. Right now it is on my
#Figure 5 from the essay

import numpy as np
import matplotlib.pyplot as plt     #to plot the graphs
plt.rcParams["figure.dpi"] = 1000   #sets a high resolution


def f(y, t): #our derivative function
    return t - y


def rk1(y0, t0, h): #RK1 method
    return y0 + h * f(y0, t0)

def rk2(y0, t0, h): #RK2 method
    k1 = h * f(y0, t0)
    k2 = h * f(y0 + k1, t0 + h)
    return y0 + 1 / 2 * (k1 + k2)
        
def rk4(y0, t0, h): #RK4 method
    k1 = h * f(y0, t0)
    k2 = h * f(y0 + k1 / 2, t0 + h / 2)
    k3 = h * f(y0 + k2 / 2, t0 + h / 2)
    k4 = h * f(y0 + k3, t0 + h)
    return y0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6


yn1 = 1  #sets y_0
yn2 = 1
yn4 = 1
rk1vals = []
rk2vals = []
rk4vals = []
finalt = 2  #sets range of time. This is needed as sometimes I wanted to show different numbers of steps


tvals1 = []
tvals2 = []
tvals4 = []
h1 = 0.1  #h value for RK1

t = 0
N1 = int(np.ceil(finalt/h1)) + 1

for i in range(N1):   #gets the RK1 values from the function
    rk1vals.append(yn1)
    yn1 = rk1(yn1, t, h1)
    tvals1.append(t)
    t += h1

t = 0
h2 = 0.2 #h value for RK2
N2= int(np.ceil(finalt/h2)) + 1

for i in range(N2):   #gets the RK2 values from the function
    rk2vals.append(yn2)
    yn2 = rk2(yn2, t, h2)
    tvals2.append(t)
    t += h2
    
t = 0
h4 = 0.4 #h value for RK4
N4 = int(np.ceil(finalt/h4)) + 1

for i in range(N4):   #gets the RK4 values from the function
    rk4vals.append(yn4)
    yn4 = rk4(yn4, t, h4)
    tvals4.append(t)
    t += h4
  
    
  
fig,axs = plt.subplots(1) #I sometimes used 2 subplots
axs.vlines(0, -10, 90, linestyle = "--", color = "purple", zorder = 0) #vertical lines to clearly show the error for 1 step (not always used)
axs.vlines(0.4, -10, 90, linestyle = "--", color = "purple", zorder = 0)


acctvals = np.linspace(-1,finalt,1000) #our t values for the "true solution" function
yvals = 2 * np.exp(-acctvals) + acctvals - 1 #calculates the true solution as we know the analytic solution to our equation
print(yvals[-1], rk1vals[-1],rk2vals[-1],rk4vals[-1]) #useful to compare the methods to the true solution
axs.plot(acctvals, yvals, label = "True Solution", linewidth = 1.5, color = "black", zorder = 1) #plots true solution
plt.scatter(0.4, 2 * np.exp(-0.4) - 0.6, color = "green", zorder = 2) #plots a scatter point of to compare the methods to this point


axs.plot(acctvals, -1 * acctvals + 1, color = "red", linestyle = "--", linewidth = 2.5, label = "Tangent at t = 0") #tangents calculated for RK2
axs.plot(acctvals, -0.2 * acctvals + 1, color = "c",linestyle = "--", linewidth = 2.5, label = "Gradient of tangent at t = 0.4, y = 0.6")


#axs.plot(tvals1, rk1vals, "--", label = "RK1 with h = " + str(h1), linewidth = 2.5, zorder = 1) #plots RK1
#plt.scatter(0.4,rk1vals[-1], color = "pink", zorder = 2) #plots RK1 scatter point
#axs.plot(tvals2, rk2vals, "--", label = "RK2 with h = " + str(h2), linewidth = 2.5, zorder = 1) #plots RK2
#plt.scatter(0.4,rk2vals[-1], color = "pink", zorder = 2) #plots RK2 scatter point
#axs.plot(tvals4, rk4vals, "--", label = "RK4 with h = " + str(h4), linewidth = 2.5, zorder = 1) #plots RK4
#plt.scatter(0.4,rk4vals[-1], color = "pink", zorder = 2) #plots RK4 scatter point
#axs.plot(tvals, error(rk1vals,tvals), label = "RK1 error with h = 1") #plots RK1 absolute error
#axs.plot(tvals, error(rk2vals,tvals), label = "RK2 error with h = 1") #plots RK2 absolute error
#axs.plot(tvals, error(rk4vals,tvals), label = "RK4 error with h = 1") #plots RK4 absolute error

axs.set_xlim(-0.05, 0.45) #sets limits and labels (not always used)
axs.set_ylim(0.58, 1.03)
axs.set_xlabel("t")
axs.set_ylabel("y")
axs.legend()

