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
u(t-dt,x) = - \alpha \cdot u(t,x-dx) + (1+2\alpha) \cdot u(t,x) - \alpha \cdot u(t,x+dx),
\end{equation*}
where again $\displaystyle \alpha \coloneqq \frac{dt}{dx^2}$.

This approach of discretization is called **<span style = "color: blue"> backward Euler method</span>.** Note that in the equation above the values $u(t,x-dx), u(t,x)$ and $u(t,x+dx)$ are unknown to us, so using the backward Euler method is not so straight-forward since it boils to solving a system of linear equations.

Although using the backward Euler method is not as direct as using the forward method, it is numerical more stable between the two so we will use the backward method for the discretization.

## Discretization of the curve shortening flow

## Implementation in Python
