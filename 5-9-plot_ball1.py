from matplotlib.pyplot import plot, show, xlabel, ylabel
from numpy import linspace
 
#freefall
vinitial = 10.
g = 9.81
t = linspace(0, (2 * vinitial)/g, 100)
 
def freefall(time):
 y = vinitial * time - 0.5 * g * time**2
 return y
 
 
xlist = t.tolist()
ylist = [freefall(i) for i in xlist]
 
plot(xlist, ylist)
xlabel('t')
ylabel('y')
show()