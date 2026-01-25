#----- REQUIRED PACKAGES -----
import numpy as np
import scipy as sp
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from numpy import pi as pi
from numpy import sin as sin
from numpy import cos as cos


#------ INITIALIZE THE CURVE -----

N = 5000 #number of points used to display the curve
dt = 0.01 #time-step
u = np.linspace(0,2*pi,N+1) #space parameters

#initial_curve is just used for the initial display
initial_curve = np.zeros((N+1,2))
x_coordinates = [17/31*sin(235/57 - 32*t) + 19/17*sin(192/55 - 30*t) + 47/32*sin(69/25 - 29*t) + 35/26*sin(75/34 - 27*t) + 6/31*sin(23/10 - 26*t) + 35/43*sin(10/33 - 25*t) + 126/43*sin(421/158 - 24*t) + 
                143/57*sin(35/22 - 22*t) + 106/27*sin(84/29 - 21*t) + 88/25*sin(23/27 - 20*t) + 74/27*sin(53/22 - 19*t) + 44/53*sin(117/25 - 18*t) + 126/25*sin(88/49 - 17*t) + 79/11*sin(43/26 - 16*t) + 
                43/12*sin(41/17 - 15*t) + 47/27*sin(244/81 - 14*t) + 8/5*sin(79/19 - 13*t) + 373/46*sin(109/38 - 12*t) + 1200/31*sin(133/74 - 11*t) + 67/24*sin(157/61 - 10*t) + 583/28*sin(13/8 - 8*t) + 
                772/35*sin(59/16 - 7*t) + 3705/46*sin(117/50 - 6*t) + 862/13*sin(19/8 - 5*t) + 6555/34*sin(157/78 - 3*t) + 6949/13*sin(83/27 -t) - 6805/54*sin(2*t + 1/145) - 
                5207/37*sin(4*t + 49/74) - 1811/58*sin(9*t + 55/43) - 63/20*sin(23*t + 2/23) - 266/177*sin(28*t + 13/18) - 2/21*sin(31*t + 7/16) for t in u]

y_coordinates = [70/37*sin(65/32 - 32*t) + 11/12*sin(98/41 - 31*t) + 26/29*sin(35/12 - 30*t) + 54/41*sin(18/7 - 29*t) + 177/71*sin(51/19 - 27*t) + 59/34*sin(125/33 - 26*t) + 49/29*sin(18/11 - 25*t) + 
                151/75*sin(59/22 - 24*t) + 52/9*sin(118/45 - 22*t) + 52/33*sin(133/52 - 21*t) + 37/45*sin(61/14 - 20*t) + 143/46*sin(144/41 - 19*t) + 254/47*sin(19/52 - 18*t) + 246/35*sin(92/25 - 17*t) + 
                722/111*sin(176/67 - 16*t) + 136/23*sin(3/19 - 15*t) + 273/25*sin(32/21 - 13*t) + 229/33*sin(117/28 - 12*t) + 19/4*sin(43/11 - 11*t) + 135/8*sin(23/10 - 10*t) + 205/6*sin(33/23 - 8*t) + 
                679/45*sin(55/12 - 7*t) + 101/8*sin(11/12 - 6*t) + 2760/59*sin(40/11 - 5*t) + 1207/18*sin(21/23 - 4*t) + 8566/27*sin(39/28 - 3*t) + 12334/29*sin(47/37 - 2*t) + 15410/39*sin(185/41 - t) - 
                596/17*sin(9*t + 3/26) - 247/28*sin(14*t + 25/21) - 458/131*sin(23*t + 21/37) - 41/36*sin(28*t + 7/8) for t in u]

x_coordinates = [x/150 for x in x_coordinates]
y_coordinates = [y/150 for y in y_coordinates]

initial_curve[:,0] = x_coordinates #initialize the x-coordinates of the initial curve
initial_curve[:,1] = y_coordinates #initialize the y-coordinates of the initial curve

#gamma_0 is being changed during the animation
gamma_0 = np.array(initial_curve[:N, :] ) #a copy of coordinates of the initial curve

time = 0.0 #set the time

#----- DISCRETIZATION OF THE CURVE SHORTENING FLOW

def csf(): #this function returns the matrix of the system we get in the discretization
    global dt, gamma_0, N 
    dGamma = np.array(gamma_0) #the denominator of the discretization of the flow
    dGamma[1:,:] = gamma_0[1:,:] - gamma_0[:-1,:] #gamma_(i+1) - gamma_i on all the points except the first one
    dGamma[0,:] = gamma_0[0,:] - gamma_0[-1,:] #gamma_(i+1) - gamma_i on the first point as well
    dGamma_squared = (dGamma*dGamma).sum(1) #|gamma_(i+1)-gamma_i|^2 for all i = 0, 1, ..., N
    temp = dGamma_squared[0]
    dGamma_squared[:-1] += dGamma_squared[1:] #|gamma_(i+1)-gamma_i|^2 + |gamma_i - gamma_(i-1)|^2 on all the points except the last
    dGamma_squared[-1] += temp #|gamma_(i+1)-gamma_i|^2 + |gamma_i - gamma_(i-1)|^2 on the last point as well
    dGamma_squared /= 2 
    min = np.amin(dGamma_squared)
    if min < 0.000001**2: #if the distance between the points becomes very small, stop the flow
        dt = 0 #stops the flow
        a = np.ones_like(dGamma_squared)
        return scipy.sparse.diags(a)
    a = dt/dGamma_squared # a = [a_0 a_1 ... a_N]
    E = sp.sparse.diags(1+2*a) #matrix which at the position (i,i) has 1+2a_i
    E += sp.sparse.diags(-a[:-1], 1,format='csc') #matrix that on the diagonal offset by 1 to the right has all the elements from -a but the last one
    E += sp.sparse.diags(-a[1:], -1,format='csc') #matrix that on the diagonal offest by 1 to the left has all the elements from -a but the first one
    E += sp.sparse.diags([-a[-1]], 1-N,format='csc') #matrix that on the bottom (N-1)st diagonal has the last element from the vector -a
    E += sp.sparse.diags([-a[0]], N-1,format='csc') #matrix that on the upper (N-1)st diagonal has the first element from the vector -a
    return E

def init(): # variable 'coordinates' represent the coordinates of the flow at time t, and 't' represents the elapsed time 
    return coordinates, t

#----- ANIMATION OF THE FLOW

def frame_generator(): #this function generates 4 extra seconds worth of frames when the flow stops
    extra_frames = 40
    extra_left = None
    global dt

    while True:
        yield None
        
        if dt == 0:
            if extra_left is None:
                extra_left = extra_frames
            else:
                extra_left -= 1

                if extra_left <= 0:
                    return


def animate(frame): #this function is actually used to animate the flow by generating frames
    global time, dt, gamma_0, initial_curve, N 
    A = csf() 
    gamma_0 = sp.sparse.linalg.spsolve(A, gamma_0) 
    initial_curve[:N, :] = gamma_0 
    initial_curve[N, :] = gamma_0[0, :] 
    time += dt 
    coordinates.set_data(initial_curve[:, 0], initial_curve[:, 1])
    if dt > 0:
        table[(0, 1)].get_text().set_text("%8.5f" % (time)) #update the time displayed in the table
    else:
         table[(0,0)].get_text().set_text("Maximal time") #update the text in the table when the flow ends
    t.set_fontsize(16)
    return coordinates, t


#----- DISPLAY OF THE COORDINATE SYSTEM

fig = plt.figure() #the figure in which everything will be displayed
plt.plot(initial_curve[:, 0], initial_curve[:, 1], "darkblue", linewidth = 0.75) #this can be commented-out if we don't want to keep the initial curve displayed during the whole animation
coordinates, = plt.plot([], [], 'r', linewidth = 1.0) #this is how we plot the curve that will be evolving
t = plt.text(2, -2.5, "")

#display of the coordinate axes
plt.xlim(-7, 7) 
plt.ylim(-7, 7)
ax = plt.gca()
ax.set_facecolor("#e6f2ff") #the background color of the plane
ax.spines["left"].set_position("zero") 
ax.spines["bottom"].set_position("zero")
ax.spines["right"].set_color("none")
ax.spines["top"].set_color("none")
ax.set_aspect('equal', adjustable='box')
plt.grid(True)

#----- DISPLAY OF THE TIME TABLE -----
table = ax.table(
    cellText=[
        ["Elapsed time", "0.00"]
    ],
    cellLoc="center",
    bbox=[0.65, 1, 0.5, 0.08], #arguments: first two are alignment args (expressed in percentages), the last two are dimensions of the frame
)

table.auto_set_font_size(False)
table.set_fontsize(8)
table.scale(1, 2)

for cell in table.get_celld().values():
    cell.set_facecolor("darkblue")
    cell.get_text().set_color("white")
    cell.get_text().set_weight("bold")

#----- ANIMATION OF THE FLOW -----
ani = FuncAnimation(fig, animate, frames=frame_generator(), init_func=init, interval=100, blit=True)
ani.save("csf-pi.mp4")

plt.show()


