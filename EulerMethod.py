###Code for DE solver by euler method, increment by small step h
import numpy as np
import matplotlib.pyplot as plt

##define derivative
def dudt(u, v, a1, b):
	return a1/(1 + v**b) - u

def dvdt(u, v, a2, g):
	return a2/(1 + u**g) - v

##define initial conditions
y0 = [0, 0]
u_0 = y0[0]
v_0 = y0[1]
a1 = 2
a2 = 4
b = 2
g = 2

tspan = [0, 12]
h = 0.1

t_0 = tspan[0]
t_f = tspan[1]
t1 = np.arange(t_0, t_f, h)
nStep = t_f/h

###Calculating
u = list()
v = list()
u.insert(0, y0[0])
v.insert(0, y0[1])
for i in range(1,int(nStep)):
	u.append(u_0 + h*dudt(u_0, v_0, a1, b))
	v.append(v_0 + h*dvdt(u_0, v_0, a2, g))
	u_0 = u[i]
	v_0 = v[i]


##Plot results
plt.scatter(t1, u, c='r')
plt.scatter(t1, v, c='b')
plt.show()


####ODE45 implementatiopn in python
from scipy.integrate import ode
import numpy as np
import matplotlib.pyplot as plt

def dudt(u, v, a1, b):
	return a1/(1 + v**b) - u

y0 = [0, 0]
t_0 = 0

r = ode(dudt).set_integrator('dopri5')
r.set_initial_value(y0[0], t_0).set_f_params(a1=2, b=2)
t1 = 12
dt = 0.1
while r.successful() and r.t < t1:
	r.integrate(r.t+dt)
	print("%g %g" % (r.t, r.y))











