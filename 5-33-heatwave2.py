from numpy import linspace, exp, sqrt, cos, log, pi
from matplotlib.pyplot import plot, show, subplots
from matplotlib.animation import FuncAnimation

k = 1E-6            # thermal diffusivity (in m**2/s)

A1 = 15             # amplitude of the daily temperature variations (in C)
P1 = 24 * 60 * 60.      # oscillation period of 24 h (in seconds)
omega1 = 2 * pi / P1   # angular freq of daily temp variations (in rad/s)
a1 = sqrt(omega1 / (2 * k))

A2 = 7                    # amplitude of yearly temperature variations (in C)
P2 = 24 * 60 * 60 * 365.  # oscillation period of 1 yr (in seconds)
omega2 = 2 * pi / P2      # angular freq of yearly temp variations (in rad/s)
a2 = sqrt(omega2 / (2 * k))

dt = P2 / 20            # time lag: 0.1 yr
t0=0
tmax = 3 * P2               # 3 year simulation
T0 = 10                 # mean surface temperature in Celsius
D = -(1 / a1) * log(0.001)  # max depth
n = 501                 # no of points in the z direction

number_of_points=int(P2/1000)

fig, ax=subplots()
timedata, temperaturedata1, temperaturedata2=[],[],[]
curve1, = plot([],[], 'ro')  
curve2, = plot([],[], 'b-')


def update(frame):
  tim=frame
  timedata.append(tim)
  temperaturedata1.append(T(1,tim))
  temperaturedata2.append(T(2,tim))
  curve1.set_data(timedata, temperaturedata1)
  curve2.set_data(timedata, temperaturedata2)
  return curve1, curve2


def init():
    ax.set_xlim(0,tmax)
    ax.set_ylim(4, 16)
    return curve1, curve2



def T(z, t):
    # T0, A, k, and omega are global variables
    return T0 + A1 * exp(-a1 * z) * cos(omega1 * t - a1 * z) + \
        A2 * exp(-a2 * z) * cos(omega2 * t - a2 * z)
print (T0-A1*exp(-a1)-A2*exp(-a2))

frames=linspace(t0, tmax, number_of_points)
print(max(T(2,frames)))
print(min(T(2,frames)))
ani=FuncAnimation(fig, update, frames=frames,init_func=init, blit=True)
show()
