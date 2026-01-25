# Numerical Simulation

This entire chapter is dedicated to Sigurd Angenent, who was kind enough to share his notes and the original code with one of the contributors of this project years ago. Most of the contents of this chapter are directly based on those materials, for which we are deeply thankful.

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
