# ESP_Cart_pendulum
mechanical system of a cart and pendulum on a spring

import numpy as np                            # libraries used
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation


%matplotlib inline


M = 5     # Mass of the cart
m = 0.1   # mass of the pendulum
k = 25    # spring constant
L = 0.5   # lenght of string
g = 9.8   # gravitational acceleration 
u1 = 0.1  # friction coefficient pendulum
u2 = 0.1   # friction coefficient cart


def der_state(t, state):                                          # Here, we are writing first order DE.
        """compute the derivative of the given state"""
        
        der = np.zeros_like(state)
        der[0] = state[1]           #    state[0] = x
        der[2] = state[3]           #    state[2] = theta
        der[1] = ((((L*m)*state[3]**2)*np.sin(state[2])) + m*g*np.sin(state[2])*np.cos(state[2]) + (u1/L)*state[3]*np.cos(state[2]) - k*state[0] - u2*state[1])/(M + m*np.sin(state[2])*np.sin(state[2]))
        der[3] = (-(m + M)*g*np.sin(state[2]) - (L*m)*(state[3]**2)*np.sin(state[2])*np.cos(state[2]) - (1 + (M/m))*(u1/L)*state[3] + np.cos(state[2])*(k*state[0] + u2*state[1]))/(L*M + L*m*np.sin(state[2])*np.sin(state[2]))
        return der
        


tf = 25 # 25 seconds of simulation
n = 1000 # number of evaluation points
dt = tf/n
T = np.linspace(0.0, tf, n+1)


state0 = ([1, 0, 1, 0]) # [x],[x'],[theta],[theta'] initial conditions

sol = integrate.solve_ivp(der_state, (0, tf), state0, t_eval=T)     #this solves both previous DE.
ang_pendulum_pos = sol.y[2]
x_quisite_position_cart = sol.y[0]


# Cartesian coordinates
x2 = x_quisite_position_cart
y2 = 0 # no motion in y direction.

x = x2 + L*np.sin(ang_pendulum_pos)
y = -L*np.cos(ang_pendulum_pos) #the y-axis points down


from matplotlib import rc
rc('animation', html='html5')


fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2.1), ylim=(-1, 1))       # setting the graph axis
ax.grid()



spring, = ax.plot([],[], linestyle = ':', lw=2)     
line_of_wall, = ax.plot([], [], 'o-', lw=7, color='red')
line_of_Table, = ax.plot([], [], 'o-', lw=7)
cart, = ax.plot([], [], 's', lw=1, markersize=15)
time_template = 'time = '
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
line_of_pendulum, = ax.plot([], [], '-o', lw=2, markersize=0)
pendulum_mass, = ax.plot([], [], 'o', lw=0, markersize=8, color='black')

def init():
    line_of_pendulum.set_data([], [])
    spring.set_data([],[])
    line_of_wall.set_data([],[])
    line_of_Table.set_data([],[])
    cart.set_data([],[])
    pendulum_mass.set_data ([],[])
    time_text.set_text('')
    return line_of_pendulum, time_text

def animate(i):
    thisLine_of_tablex = [-2, 2]                # x coordinates of table are fixed
    thisLine_of_tabley = [-0.035,-0.035,]       # y coordinates of table are fixed
    thisSpringx = [x2[i], 2]                 # x coordinates of spring, one to cart one to wall
    thisSpringy = [0.05, 0.05]               # y coordinares of spring, both fixed
    thisLine_of_wallx = [2, 2]                  # x coordinates of wall are fixed
    thisLine_of_wally = [-0.035, 1]             # y coordinates of wall are fixed
    thisLine_of_pendulumx = [x2[i], x[i]]    # x coordinates of pendulum line 'L' are moving, one with the cart one with the pendulum mass
    thisLine_of_pendulumy = [0.035, y[i]]    # one y coordinate of pendulum line 'L' is fixed with the cart, one is moving with the pendulum mass 
    thisCartx = [x2[i],x2[i]]                   # x coordinates of the cart are moving, the cart is moving.
    thisCarty = [0.05, 0.05]                    # y coordinates of the cart are fixed
    thispendulumx = [x[i], x[i]]             # x coordinates of the pendulum mass 'm' are moving
    thispendulumy = [y[i], y[i]]             # y coordinates of the pendulum mass 'm' are moving
  

    line_of_pendulum.set_data(thisLine_of_pendulumx, thisLine_of_pendulumy)
    spring.set_data(thisSpringx, thisSpringy)
    line_of_wall.set_data(thisLine_of_wallx, thisLine_of_wally)
    line_of_Table.set_data(thisLine_of_tablex, thisLine_of_tabley)
    cart.set_data(thisCartx, thisCarty)
    pendulum_mass.set_data(thispendulumx, thispendulumy)
    

    time_text.set_text(time_template + '{:4.1f}'.format(i*dt) + 's')
    return line_of_pendulum, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(T)),
                              interval=20, blit=True, init_func=init)
plt.close(fig)


ani

