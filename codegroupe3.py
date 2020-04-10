from math import pi, cos
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D     # useful for the phase diagram

# Parameters:

# sine function of our air temperature :
w = 2 * pi / 6      # pulsation
phi = 3             # phase, in radiant (to start when the temperature is rising)
Am = 20             # Amplitude of Air temperature

kAI = 0.001         # Impact of air temperature on ice
kAS = 0.003         # Impact of air temperature on snow
kAG = 0.6           # Impact of air temperature on ground temperature

kB = 0.05           # Impact of Bisons

kGS = 0.003         # Impact of ground temperature on snow, lower than the one below
kGI = 0.004         # Impact of ground temperature on ice

kIG = 0.05          # Impact of ice on ground temperature

# we assumed that the ground temperature try to go through air temperature
B = 20               # density of bizon


ListA = []
Listt = []


def fonction(y, t, args):
   # input variable, we added 264 to adjust from when the temperature starts as we
   # wanted 284K fo the maximum and 234 K for the minimum
    A = (Am) * cos(w * t + phi) + 264
    ListA.append(A)
    Listt.append(t)

    S, I, G = y
    kAS, kB, kGS, kAI, kGI, kAG = args

    # differential equation :
    dSdt = -kAS * A * S - kGS * G * S - kB * B * S + kAS * S
    dIdt = -kAI * A * I - kGI * G * S + kB * B * S
    dGdt = kAG * (A - G) - kIG * I

    dydt = [dSdt, dIdt, dGdt]
    return dydt


y0 = [400, 300, 274]        # with snow quantity = 400 cm, ice quantity = 300 cm, air temperature = 274 K
t = np.linspace(1, 12, 48)
sol = odeint(fonction, y0, t, args=([kAS, kB, kGS, kAI, kGI, kAG],))


# plot :
fig = plt.figure()
fig.subplots_adjust(hspace=0.7, wspace=0.7)
fig.add_subplot(2, 1, 1)
plt.plot(Listt, ListA, label="Air ")
plt.plot(t, sol[:, 2], "g", label="Ground Temperature")
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.45), ncol=3, fancybox=True, shadow=True)
plt.xlabel("Time in months")
plt.ylabel("Air temperature in Kelvin")
plt.title("Air Temperature and Ground Temperature variation")

fig.add_subplot(2, 1, 2)
plt.plot(t, sol[:, 0], "r", label="Snow quantity")
plt.plot(t, sol[:, 1], "b", label="Ice quantity")
plt.legend(loc='best')
plt.xlabel("Time in month")
plt.ylabel("Snow and Ice layer thickness \n in cm")
plt.title("Snow and Ice quantity variation")

plt.show()


# phase diagram
figure = plt.figure()
ax = figure.gca(projection='3d')
ax.plot(sol[:, 0], sol[:, 1], sol[:, 2])
ax.set_xlabel('dSdt\nVariation of snow quantity in time')
ax.set_ylabel('dIdt\nVariation of ice quantity in time')
ax.set_zlabel('dGdt\nVariation of ground temperature in time')
plt.show()
