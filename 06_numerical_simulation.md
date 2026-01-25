# Numerical Simulation

*This entire chapter is dedicated to Sigurd Angenent, who was kind enough to share his notes and the original code with one of the contributors of this project years ago. Most of the contents of this chapter are directly based on those materials, for which we are deeply thankful.*

## Discretization of the 1-dimensional heat equation
Consider $1$-dimensional heat equation given by
\begin{equation*}
u_t = u_{xx},
\end{equation*}
where $t$ represents time and $x$ represents the space parameter.

In order to discretize the heat equation, we will consider the values of the solution $u(t,x)$ on a grid
\begin{equation*}
G = \left \{(M \cdot dt, N \cdot dx) : M,N \in \mathbb{N} \right \},
\end{equation*}
where $dt$ and $dx$ represent the distance between consecutive time and space parameters on the grid $G$. 

### Forward Euler Method

The derivative $u_t$ can be discretized as
\begin{equation*}
u_t(t,x) = \frac{u(t+dt, x) - u(t,x)}{dt}.
\end{equation*}

As for the space derivative, the Taylor expansion gives us the following:
\begin{equation*}
\begin{split}
u(t, x + dx ) &= u(t,x) + u_x(t,x) \cdot dx + \frac{1}{2}u_{xx}(t,x) \cdot dx^2 + \mathcal{O}(dx^2) \\
u(t, x-dx) &= u(t,x) -u_x(t,x) \cdot dx + \frac{1}{2}u_{xx}(t,x) \cdot dx^2 + \mathcal{O}(dx^2)
\end{split}
\end{equation*}
Adding these two equations we get
\begin{equation*}
u_{xx}(t,x) = \frac{u(t,x+dx) -2u(t,x) + u(t, x-dx)}{dx^2} + \mathcal{O}(dx^2).
\end{equation*}
Therefore, one possible discretization of the heat equation is given by
\begin{equation*}
\frac{u(t+dt, x) - u(t,x)}{dt} = \frac{u(t,x+dx) -2u(t,x) + u(t, x-dx)}{dx^2}
\end{equation*}

As we will see in the section on implementation in Python, the solution $u(t, \cdot)$ at time $t$ will be known to us and we will be using it to find the solution $u(t + dt, \cdot)$ at time $t+dt$. From the equation above we see that
\begin{equation*}
u(t+dt,x) = \alpha \cdot u(t, x-dx) + (1-2\alpha) \cdot u(t,x) + \alpha \cdot u(t,x+dx),
\end{equation*}
where $\displaystyle \alpha \coloneqq \frac{dt}{dx^2}$.

This approach of discretization is called **<span style = "color: blue"> forward Euler method</span>,** and using it is very straight-forward since computing $u(t+dt,x)$ boils down to solving a single equation.

### Backward Euler Method
A variation of the forward Euler method can be done by discretizing the time derivative $u_t$ in the following way:
\begin{equation*}
u_t(t,x) = \frac{u(t,x) - u(t-dt,x)}{dt}.
\end{equation*}
Analogous computations as before yield the following discretization of the heat equation:
\begin{equation*}
\label{BWD}
u(t-dt,x) = - \alpha \cdot u(t,x-dx) + (1+2\alpha) \cdot u(t,x) - \alpha \cdot u(t,x+dx),
\end{equation*}
where again $\displaystyle \alpha \coloneqq \frac{dt}{dx^2}$.

This approach of discretization is called **<span style = "color: blue"> backward Euler method</span>.** Note that in the equation above the values $u(t,x-dx), u(t,x)$ and $u(t,x+dx)$ are unknown to us, so using the backward Euler method is not so straight-forward since it boils to solving a system of linear equations.

Although using the backward Euler method is not as direct as using the forward method, it is numerical more stable between the two so we will use the backward method for the discretization.

### Pseudo-algorithm for the backward Euler method
Now we want to see how to explicitly use the backward Euler method in order to compute the values of the solution $u(t, \cdot)$ at time $t$ given the values $u(t-dt, \cdot)$ at time $t-dt$.

For simplicity, assume that the space paramter $x$ of the heat equation belongs to the segment $[0,1]$ and assume that our solution is $1$-periodic, i.e. that $u(t,x+1) = u(t,x)$. Consider the subdivision of the segment $[0,1]$ given by the points
\begin{equation*}
\left \{0, \frac{1}{N}, \frac{2}{N}, \dots, \frac{N-1}{N}, 1 \right \},
\end{equation*}
where $N \in \mathbb{N}$ is some natural number that corresponds to the number of points used to discretize the space parameter. [The backward Euler method](#BWD) tells us that in order to compute $u(t,\cdot)$ given $u(t-dt, \cdot)$, we need to solve the following system of linear equations:
\begin{equation*}
\begin{bmatrix} 
1+2\alpha & -\alpha & 0 & 0 & \dots & 0 & - \alpha \\
-\alpha & 1+2\alpha & -\alpha & 0 & \dots & 0 & 0 \\
0 & -\alpha & 1+2\alpha & - \alpha & \dots & 0 & 0 \\
\vdots & \vdots & \vdots & \vdots & \dots & \vdots & \vdots \\
-\alpha & 0 & 0 & 0 & \dots & -\alpha & 1+2\alpha
\end{bmatrix}
\cdot
\begin{bmatrix}
u(t,0) \\ u(t,\frac{1}{N}) \\ u(t,\frac{2}{N}) \\ \vdots \\ u(t, \frac{N-1}{N})
\end{bmatrix} =
\begin{bmatrix}
u(t-dt,0) \\ u(t-dt,\frac{1}{N}) \\ u(t-dt,\frac{2}{N}) \\ \vdots \\ u(t-dt, \frac{N-1}{N})
\end{bmatrix}
\end{equation*}
Now, let
\begin{equation*}
\begin{split}
E &= \begin{bmatrix}
2 & -1 & 0 & 0 & \dots & 0 & -1 \\
-1 & 2 & -1 & 0 & \dots & 0 & 0 \\
0 & -1 & 2 & -1 & \dots & 0 & 0 \\
\vdots & \vdots & \vdots & \vdots & \dots & \vdots & \vdots \\
-1 & 0 & 0 & 0 & \dots & -1 & 2 
\end{bmatrix} \\
A&= I + \alpha \cdot E \\
U(t) &= \begin{bmatrix} u(t,0) \\ u(t, \frac{1}{N}) \\ u(t, \frac{2}{N}) \\ \vdots \\ u(t, \frac{N-1}{N}) \end{bmatrix}
\end{split}
\end{equation*}
Using this notation, the system above can be written as
\begin{equation*}
\begin{split}
A &\cdot U(t) = U(t-dt) \\
\implies U(t) &= A^{-1} \cdot U(t-dt).
\end{split}
\end{equation*}

This gives us the following pseudo-code on how to discretize the heat equation using the backward Euler method:

:::{note} ⚙️ The backward Euler method algorithm
:icon: false
**Inputs:**
- number of space parameters $\dots N$
- time step $\dots dt$
- initial condition $\dots u(0,x)$

**Output:**

1. Compute $\displaystyle dx = \frac{1}{N}$
2. Compute the matrix $A \in M_n(\mathbb{R})$
3. Initialize the time to $t = 0$ and set $U = u(0,x)$
4. **Repeat**:

    - compute the solution at the time $t+dt$ by calculating $A^{-1} \cdot U$
    - update the time $t + dt$
:::

## Discretization of the curve shortening flow

In this section we are going to discretize the curve shortening flow. Let $\gamma_0 : I \rightarrow \mathbb{R}^2$ be a regular parametrized curve and let $\gamma(t,u)$ be its curve shortening flow, where $u$ is a general parameter. The curve shortening flow then boils down to the following equation:
\begin{equation*}
\frac{\partial \gamma}{\partial t} = \frac{1}{\lvert \gamma' \rvert} \cdot \frac{\partial}{\partial u} \left(\frac{1}{\lvert \gamma'\rvert} \frac{\partial \gamma}{\partial u}\right).
\end{equation*}
Explicitly computing the derivate, we get that the curve shortening flow equation is 
\begin{equation*}
\frac{\partial \gamma}{\partial t} = \frac{\gamma_{uu}}{\lvert \gamma_u \rvert^2} - \frac{\langle \gamma_u, \gamma_{uu} \rangle}{\lvert \gamma_u \rvert^4}\gamma_u.
\end{equation*}
Due to the invariance of the curve shortening flow under tangential perturbations, the curve shortening flow is equivalent to the equation
\begin{equation*}
\frac{\partial \gamma}{\partial t} = \frac{\gamma_{uu}}{\lvert \gamma_u \rvert^2}.
\end{equation*}
This is the equation that we are going to discretize, and we are going to do so using the [backward Euler method](#BWD).

Discretization of the time-derivative and second space-derivative are given as before:

\begin{equation*}
\begin{split}
\frac{\partial \gamma}{\partial t} &= \frac{\gamma(t,u) - \gamma(t-dt, u)}{dt} \\
\gamma_{uu}(t,u) &= \frac{\gamma(t,u+du) - 2\gamma(t,u) + \gamma(t, u-du)}{du^2}
\end{split}
\end{equation*}

In order to discretize the norm $\lvert \gamma_u \rvert^2$, we have made the following choice:

\begin{equation*}
\lvert \gamma_u (t,u) \rvert^2 = \frac{\lvert \gamma(t-dt,u+du) - \gamma(t-dt,u) \rvert^2 + \lvert \gamma(t-dt,u) - \gamma(t-dt,u-du)\rvert^2}{2du^2}.
\end{equation*}

Now the discretized equation for the curve shortening flow is given by 

\begin{equation*}
\label{CSF-discrete}
\frac{\gamma(t,u) - \gamma(t-dt,u)}{dt} = 2 \cdot \frac{\gamma(t,u+du) - 2\gamma(t,u) + \gamma(t, u-du)}{\lvert \gamma(t-dt,u+du) - \gamma(t-dt,u) \rvert^2 + \lvert \gamma(t-dt,u) - \gamma(t-dt,u-du)\rvert^2}.
\end{equation*}

If we define

\begin{equation*}
\alpha \coloneqq \frac{2\cdot dt}{\lvert \gamma(t-dt,u+du) - \gamma(t-dt,u) \rvert^2 + \lvert \gamma(t-dt,u) - \gamma(t-dt,u-du)\rvert^2},
\end{equation*}

then we get a system of linear equations that is analogous to the one for the heat equation:

\begin{equation*}
-\alpha \cdot \gamma(t,u-du) + (1+2\alpha) \cdot \gamma(t,u) - \alpha \cdot \gamma(t, u + du) = \gamma(t-dt, u).
\end{equation*}
Hence, we can use the same pseudo-algorithm to find its solution:

:::{note} ⚙️ The backward Euler method algorithm for the curve shortening flow
:icon: false
**Inputs:**
- number of space parameters $\dots N$
- time step $\dots dt$
- initial condition $\dots \gamma_0(u)$

**Output:**

1. Compute $\displaystyle du = \frac{1}{N}$
2. Compute the matrix $A \in M_n(\mathbb{R})$
3. Initialize the time to $t = 0$ and set $\gamma = \gamma(0,u)$
4. **Repeat**:

    - compute the solution at the time $t+dt$ by calculating $A^{-1} \cdot \gamma$
    - update the time $t + dt$
:::

## Implementation in Python

In this section, we implement the pseudo-algorithm for the curve shortening flow in Python and we generate the animations of the flow for various initial curves. Regarding the display, we want to achieve the following:
- display the initial curve
- display the curve shortening flow, keeping the initial curve fixed
- display the elapsed time of the flow

First of all, we need to initialize number of space parameters, the time step and the initial curve. In order to present the basic outline of the code, we are going to discretize the curve shortening flow of the circle with centre $(0,0)$ and radius equal to $2$.

```{code} python
:caption: The first part of code in which we initialize the number of points, the time step and the initial curve.
#------ INITIALIZE THE CURVE -----

N = 5000 #number of points used to display the curve
dt = 0.01 #time-step
u = np.linspace(0,1,N+1) #space parameters

#initial_curve is just used for the initial display
initial_curve = np.zeros((N+1,2)) 
initial_curve[:,0] = 2*np.cos(2*np.pi*u) #initialize the x-coordinates of the initial curve
initial_curve[:,1] = 2*np.sin(2*np.pi * u) #initialize the y-coordinates of the initial curve

#gamma_0 is being changed during the animation
gamma_0 = np.array(initial_curve[:N, :] ) #a copy of coordinates of the initial curve

time = 0.0 #set the time
```

Secondly, the following code allows us to compute the matrix $A = I + \alpha \cdot E$ of the system of linear equations we get when discretizing the curve shortening flow as described previously.

```{code} python
:caption: The part of the code in which we compute the matrix of the system of linear equations we get when discretizing the curve shortening flow.
def csf(): #this function returns the matrix of the system we get in the discretization
    global dt, gamma_0, N #variables defined globally

    #the denominator of the discretization of the flow
    dGamma = np.array(gamma_0) 
    dGamma[1:,:] = gamma_0[1:,:] - gamma_0[:-1,:] #gamma_(i+1) - gamma_i on all the points except the first one
    dGamma[0,:] = gamma_0[0,:] - gamma_0[-1,:] #gamma_(i+1) - gamma_i on the first point as well
    dGamma_squared = (dGamma*dGamma).sum(1) #|gamma_(i+1)-gamma_i|^2 for all i = 0, 1, ..., N
    temp = dGamma_squared[0]
    dGamma_squared[:-1] += dGamma_squared[1:] #|gamma_(i+1)-gamma_i|^2 + |gamma_i - gamma_(i-1)|^2 on all the points except the last
    dGamma_squared[-1] += temp #|gamma_(i+1)-gamma_i|^2 + |gamma_i - gamma_(i-1)|^2 on the last point as well
    dGamma_squared /= 2

    #if the distance between the points becomes very small, stop the flow
    min = np.amin(dGamma_squared)
    if min < 0.000001**2: 
        dt = 0 #stops the flow
        a = np.ones_like(dGamma_squared)
        return scipy.sparse.diags(a)
    
    #else, return the matrix A = I + a*E
    a = dt/dGamma_squared # a = [a_0 a_1 ... a_N]
    E = sp.sparse.diags(1+2*a) #matrix which at the position (i,i) has 1+2a_i
    E += sp.sparse.diags(-a[:-1], 1,format='csc') #matrix that on the diagonal offset by 1 to the right has all the elements from -a but the last one
    E += sp.sparse.diags(-a[1:], -1,format='csc') #matrix that on the diagonal offest by 1 to the left has all the elements from -a but the first one
    E += sp.sparse.diags([-a[-1]], 1-N,format='csc') #matrix that on the bottom (N-1)st diagonal has the last element from the vector -a
    E += sp.sparse.diags([-a[0]], N-1,format='csc') #matrix that on the upper (N-1)st diagonal has the first element from the vector -a
    return E
```

Lastly, we want to animate the flow and in order to do so, we use the *Animation* class from *Matplotlib* library. The main part here is to setup the function for generating frames, and that is done with the following code:

```{code} python
:caption: The part of the code in which we setup the function for generating frames of the animation and we display the solution of the flow at each time $t$.
#----- ANIMATION OF THE FLOW

# variable 'coordinates' represent the coordinates of the flow at time t, and 't' represents the elapsed time
def init():  
    return coordinates, t

#this function generates 4 extra seconds worth of frames when the flow stops
def frame_generator(): 
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

#this function is actually used to animate the flow by generating frames
def animate(frame): 
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
plt.xlim(-3, 3) 
plt.ylim(-3, 3)
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
ani.save("csf-circle.mp4")

plt.show()
```

Below you can find the entire code.

```{code} python
:caption: The full Python code used for the discretizaion and the display of the curve shortening flow of a circle with centre $(0,0)$ and radius $2$.
:label: csf-circle

#----- REQUIRED PACKAGES -----
import numpy as np
import scipy as sp
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation 

#------ INITIALIZE THE CURVE -----

N = 5000 #number of points used to display the curve
dt = 0.01 #time-step
u = np.linspace(0,1,N+1) #space parameters

#initial_curve is just used for the initial display
initial_curve = np.zeros((N+1,2)) 
initial_curve[:,0] = 2*np.cos(2*np.pi*u) #initialize the x-coordinates of the initial curve
initial_curve[:,1] = 2*np.sin(2*np.pi * u) #initialize the y-coordinates of the initial curve

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
plt.xlim(-3, 3) 
plt.ylim(-3, 3)
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
ani.save("csf-circle.mp4")

plt.show()
````

Running [the code above](#csf-circle) we get the following animation:

:::{figure} ./video/circle.mp4
Running the code above produces the animation of the curve shortening flow of the circle with center $(0,0)$ and radius equal to $2$.
:::

In the code above, you can use any initial curve that you'd like to produce all kinds of animations. Below we showcase two such possible examples.

:::{figure} ./video/csf-pi.mp4
The Gage-Hamilton-Graysong Theorem in full effect - under the curve shortening flow, a curve that is not initially convex becomes convex and then shrinks to a point in finite time.
:::

:::{figure} ./video/csf-ying-yang.mp4
This animation also illustrates a non-trivial fact about the curve shortening flow: if the initial curve is embedded, then the solution $\gamma(t, \cdot)$ is also embedded at all times $t$, meaning that there are no self-intersections during the evolution. In other words, the singularities can only occur when the curve shrinks to a point.
:::