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

N = 9999 #number of points used to display the curve
dt = 0.01 #time-step
u = np.linspace(-17,20,N+1) #space parameters

# set everything up for the ying yang curve
x_coordinates = np.zeros(N+1)
y_coordinates = np.zeros(N+1)

x_coordinates_neg_branch = []
y_coordinates_neg_branch = []
x_coordinates_pos_branch = []
y_coordinates_pos_branch = []

counter = 0
for x in u:
    if (x < 0 and x + 1/(N+1) < 0):
        x_coordinates_neg_branch.append(-np.sqrt(-x)*cos(-x))
        y_coordinates_neg_branch.append(-np.sqrt(-x)*sin(-x))
    if (x > 0 and counter == 0):
        x_coordinates_pos_branch.append(0)
        y_coordinates_pos_branch.append(0)
        counter = counter + 1
        x_coordinates_pos_branch.append(np.sqrt(x)*cos(x))
        y_coordinates_pos_branch.append(np.sqrt(x)*sin(x))
    if (x > 0 and counter != 0):
        x_coordinates_pos_branch.append(np.sqrt(x)*cos(x))
        y_coordinates_pos_branch.append(np.sqrt(x)*sin(x))

x_coordinates = x_coordinates_pos_branch + x_coordinates_neg_branch
y_coordinates = y_coordinates_pos_branch + y_coordinates_neg_branch

#initial_curve is just used for the initial display
initial_curve = np.zeros((N+3,2))
initial_curve[:,0] = x_coordinates #initialize the x-coordinates of the initial curve
initial_curve[:,1] = y_coordinates #initialize the y-coordinates of the initial curve

#gamma_0 is being changed during the animation
gamma_0 = np.array(initial_curve[:N+2, :] ) #a copy of coordinates of the initial curve

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
    E += sp.sparse.diags([-a[-1]], 1-(N+2),format='csc') #matrix that on the bottom (N-1)st diagonal has the last element from the vector -a
    E += sp.sparse.diags([-a[0]], (N+2)-1,format='csc') #matrix that on the upper (N-1)st diagonal has the first element from the vector -a
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
    initial_curve[:N+2, :] = gamma_0 
    initial_curve[N+2, :] = gamma_0[0, :] 
    time += dt 
    coordinates.set_data(initial_curve[:, 0], initial_curve[:, 1])
    if dt > 0:
        table[(0, 1)].get_text().set_text("%8.5f" % (time)) #update the time displayed in the table
    else:
         table[(0,0)].get_text().set_text("Maximal time") #update the text in the table when the flow ends
    t.set_fontsize(16)

    plt.cla() #clear axis
    plt.xlim(-7,7)
    plt.ylim(-7,7)
    ax = plt.gca()
    ax.set_facecolor("#e6f2ff") #the background color of the plane
    ax.spines["left"].set_position("zero") 
    ax.spines["bottom"].set_position("zero")
    ax.spines["right"].set_color("none")
    ax.spines["top"].set_color("none")
    ax.set_aspect('equal', adjustable='box')
    #plt.grid(True)
    plt.fill_between(initial_curve[:,0], initial_curve[:,1])

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
#plt.grid(True)

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
ani.save("csf-ying-yang.mp4")

plt.show()


