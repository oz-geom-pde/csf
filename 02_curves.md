# Curve Theory

# Parametrisation
## Definition of a Parametrised Curve
:::{prf:definition}
:label: parametrised-curve

A parametrised curve is a pair of functions $x(t), y(t)$ that describe a curve in the $\mathbb{R}^2$ plane as:
\begin{equation*}
c(t) = (x(t), y(t))
\end{equation*}

where $t \in \mathbb{R}$.
:::

## Definition of a Regular Parametrisation
:::{prf:definition}
:label: regular-parametrisation

A curve is called a regular curve if it is $C^1$, the velocity vector $c'(t)$ is never zero:
\begin{equation*}
c'(t) \neq 0 \quad \text{for all } t
\end{equation*}
:::

# Arc Length
## Arc Length Parameter
:::{prf:definition}
:label: arc-length-parameter

The arc length parameter $s(t)$ measures the distance travelled along a curve from a fixed starting point.

Given a regular curve $c(t)$, the arc length from $t = a$ to $t$ is:

\begin{equation*}
s(t) = \int_a^t \|c'(u)\| \, du
\end{equation*}
:::

## Arc Length Parametrisation
:::{prf:definition}
:label: arc-length-parametrisation

Let
\begin{equation*}
s: u \in [a, b] \rightarrow [c, d]
\end{equation*}

and its inverse
\begin{equation*}
s^{-1}: v \in [c, d] \rightarrow [a, b]
\end{equation*}

Then we can write
\begin{equation*}
u = s^{-1}(v)
\end{equation*}

Thus we can define the curve $c$ as $\tilde{c}$ or $\hat{c}$
\begin{equation*} \tilde{c}(v) = c(s^{-1}(v)) \end{equation*}
\begin{equation*} \|\partial_v \tilde{c}\| = \|\partial_v \hat{c}(s^{-1}(v))\| \end{equation*}
\begin{equation*} = \|\partial_v s^{-1} \cdot \partial_u \hat{c}\| \quad \text{(chain rule)} \end{equation*}
\begin{equation*} = \left\| \frac{1}{\partial_u s} \cdot \partial_u \hat{c} \right\| \end{equation*}
\begin{equation*} = \left\| \frac{\partial_u c}{\|\partial_u c\|} \right\| \end{equation*}
\begin{equation*} = 1 \end{equation*}


In terms of $v$:
\begin{equation*}
s(v) = \int_c^v \|\partial_w \tilde{c}(w)\|\, dw = \int_c^v dw = v - c
\end{equation*}

Arc length between $v_1$ and $v_2$:
\begin{equation*}
\ell(v_1, v_2) = \int_{v_1}^{v_2} \|\partial_w \tilde{c}(w)\|\, dw = v_2 - v_1
\end{equation*}
:::

## Arc Length of a Segment
:::{prf:definition}
:label: arc-length-segment

The arc length of a segment of a curve from $t = t_1$ to $t = t_2$ where $a,b \in \mathbb{R}$ is given by:

\begin{equation*}
\ell(t_1,t_2) = \int_{t_1}^{t_2} \|c'(t)\| \, dt
\end{equation*}
:::

## Total Arc Length
:::{prf:definition}
:label: total-arc-length

The total arc length of a curve $c(t)$ over the full interval of $c(t): t \in [a, b]$ is:

\begin{equation*}
L = \int_a^b \|c'(t)\| \, dt
\end{equation*}

This gives the total distance travelled along the curve from start to end.
:::

# Area
## Green's Theorem
:::{prf:lemma}
:label: greens-theorem

Let $D \subseteq \mathbb{R}^2$ be an open region with boundary $C$, Then by Green’s Theorem:

\begin{equation*}
\oint_C P \,dx + Q \, dy= \iint_u \partial_x Q - \partial_y P \, dA,
\end{equation*}
where $C$ is oriented positively (counter-clockwise) so $D$ is on the left when travelling around $C$.
:::

## Definition of an Enclosed Area
:::{prf:definition}
:label: enclosed-area
The **area enclosed** by $C$ is given by:
\begin{equation*}
A = \iint_u \, dA
\end{equation*}

We can find $P, Q$ such that $\partial_x Q - \partial_y P = 1$
\begin{equation*} \text{Let } \partial_x Q = \frac{1}{2} \rightarrow Q = \frac{x}{2} \end{equation*}
\begin{equation*} \text{Let } \partial_y P = -\frac{1}{2} \rightarrow P = -\frac{y}{2} \end{equation*}

\begin{equation*}
\begin{split}
\text{Area } &= \iint_u \, dA \\
&= \iint_u \partial_x Q - \partial_y P \, dA \\
& = \oint_C P \,dx + Q \, dy \quad \text{by Green's Theorem} \\
& = \oint_C -\frac{y}{2} \,dx + \frac{x}{2} \, dy \\
\end{split}
\end{equation*}
where $c(u) = (x(u),y(u)).$
\begin{equation*}
\begin{split}
A &= \frac{1}{2} \int_a^b -yx' + xy'\,du \\
&= \frac{1}{2} \int \langle (x,y), (y', -x') \rangle \, du
\end{split}
\end{equation*}

For simplicity, and application to curvature:
\begin{equation*} \text{Let } v=\|c'\|, \end{equation*}
\begin{equation*} T = \frac{c'}{\|c'\|} = \frac{(x', y')}{v} \end{equation*}
\begin{equation*} N = R_{\pi/2}(T) = \frac{(-y', x')}{v} \, du \end{equation*}

Then:
\begin{equation*} (y', -x') = vN \end{equation*}
Therefore, continuing the derivation of area:
\begin{equation*}
\begin{split}
A &= \frac{1}{2} \int \langle c, N \rangle v \, du \\
&= \frac{1}{2} \int \langle c, N \rangle\, ds
\end{split}
\end{equation*}

Therefore:
\begin{equation*}
\partial_t A = \frac{1}{2} \int \langle \partial_t c, N \rangle
+ \langle c, \partial_t N \rangle \, ds
+ \frac{1}{2} \int \langle c, N \rangle \, ds
\end{equation*}
:::



# Simple, Closed Curves
## Definition of Simple, Closed Curves
:::{prf:definition}
:label: simple-closed-curve

A curve $c: [a, b] \in \mathbb{R}^2$ is called a simple closed curve if:

- It is smooth on $[a, b]$: $c$ is twice continuously differentiable on $[a, b]$
- It is injective on $(a, b)$: $c(t_1) \ne c(t_2)$ for all $a < t_1 < t_2 < b$
- The endpoints match smoothly:
\begin{equation*} 
\begin{split}
c(a) &= c(b) \\
c'(a) &= c'(b) \\
c''(a) &= c''(b)
\end{split}
\end{equation*}
This means the curve traces out a loop without crossing itself, and joins up smoothly at the endpoints.
:::

# Curvature
## Definition of Curvature
:::{prf:definition}
:label: curvature-definition

Let $c(t) = (x(t), y(t))$ be a regular smooth curve parametrised by arc-length $s$. The **curvature** $\kappa(s)$ at a point on the curve is defined as:
\begin{equation*}
\kappa(s) = \langle \partial_s T, N \rangle
\end{equation*}

Here:
- $T = \partial_s c$ is the unit tangent vector.
- $N = R_{\pi/2} (T)$ is the unit normal vector.

This gives $\kappa$ a signed value, depending on whether the curve bends to the left (positive) or right (negative).
:::

:::{prf:remark}
:label: curvature-interpretation

Curvature measures the rate at which the tangent vector turns as a particle moves along the curve. It captures how intense the "bend" of the curve is at each point:
- If $\kappa(s) = 0$, the curve is straight.
- As $\|\kappa(s)\|$ increases, a sharper curve is indicated at that particular point.
:::

## Dependence on Unit Normal

:::{prf:remark}
:label: unit-normal-dependence

The curvature vector $\partial_s T$ points in the direction of the unit normal vector $N$, i.e.:

\begin{equation*}
\partial_s T = \kappa N
\end{equation*}

However, since $N$ is a general case, and broadly defined, both $N$ and $-N$ are orthogonal to $T$, this introduces a **sign ambiguity** in the curvature.

To resolve this, we've chosen the orientation of the plane so that counterclockwise is positive. From here, we defined $N = R_{\pi/2}(T)$ so that the pair $(T, N)$ is **positively oriented** (i.e. turning from $T$ to $N$ is counter-clockwise). Therefore $\kappa$ is well-defined with a consistent sign/direction.
:::



# Frenet-Serret Formulae
## Product Rule
:::{prf:lemma}
:label: Product rule

\begin{equation*} \text{If } X, Y: (a,b) \rightarrow \mathbb{R}^2 \end{equation*}
Where $X(t) = (x_1(t), x_2(t))$ and $Y(t) = (y_1(t), y_2(t))$
\begin{equation*}
\text{Then } \frac{d}{dt} \langle X(t), Y(t) \rangle =
\langle \frac{d}{dt} X(t), Y(t) \rangle + \langle X(t), \frac{d}{dt} Y(t) \rangle
\end{equation*}

\begin{equation*}
\text{Using dot product :}
\quad \frac{d}{dt}(X \cdot Y) = \frac{d}{dt}X \cdot Y + X \cdot \frac{d}{dt}Y
\end{equation*}

\begin{equation*}
(X \cdot Y)' = X' \cdot Y + X \cdot Y'
\end{equation*}
:::

:::{prf:proof}
:label: Product rule
\begin{equation*}
\begin{split}
\left( X \cdot Y\right)' &= (x_1 y_1 + x_2 y_2)' \\
&= (x_1 y_1)' + (x_2 y_2)' \\
&= x_1' y_1 + x_1 y_1' + x_2' y_2 + x_2 y_2' \\
&= (x_1', x_2') \cdot (y_1, y_2) + (x_1, x_2) \cdot (y_1', y_2') \\
&= X' \cdot Y + X \cdot Y'. \quad \blacksquare
\end{split}
\end{equation*}
:::

## Orthogonality of the derivative of a constant length vector field

:::{prf:lemma}
:label:
If $\| V \| \equiv \alpha,$ then $\partial_s V \perp V.$ In particular,
\begin{equation*}
\partial_s T \perp T, \quad \partial_s N \perp N.
\end{equation*}
:::

:::{prf:proof}
:label:

Since $\|V \| \equiv \alpha$ (i.e. for every $u$, $\|{V(u)}\| = \alpha$), we have
\begin{equation*}
\begin{split}
0 &= \partial_u \left( \alpha^2 \right) \\
&= \partial_u \|V(u)\|^2 \\
&= \partial_u (V \cdot V)' \quad \text{ here, } V' \text{ refers to } \partial_u V\\
&= (V' \cdot V + V \cdot V') \quad \text{by product rule}\\
&= 2(V' \cdot V).
\end{split}
\end{equation*}
Then we must have $V' \cdot V = 0.$ In particular, for $V = T$ and $V = N$ we have
\begin{equation*}
T' \cdot T = 0, \quad N' \cdot N = 0.
\end{equation*}
Taking the derivative with respect to the arc-length parameter yields:
\begin{equation*}
\therefore \partial_s T \cdot T = \left(\frac{1}{v} T'\right) \cdot T = \frac{1}{v} \left( T' \cdot T\right) = 0 \quad \blacksquare
\end{equation*}
:::

## Frenet-Seret Equations
:::{prf:theorem}
:label: frenet-serret-theorem

\begin{equation*}
\begin{split}
\partial_s T &= \kappa N \\
\partial_s N &= - \kappa T.
\end{split}
\end{equation*}
Or, in matrix form:
\begin{equation*}
\partial_s \begin{pmatrix} T \\ N \end{pmatrix}= \begin{pmatrix} 0 & \kappa \\ -\kappa & 0\end{pmatrix} \begin{pmatrix} T \\ N \end{pmatrix}
\end{equation*}
:::


:::{prf:proof}
:label: frenet-serret-proof
1. $ \partial_s T = \kappa N $
\begin{equation*}\text{By definition: } \quad \kappa = \langle \partial_s T, N \rangle \end{equation*}
\begin{equation*} \text{Since } \partial_s T \perp T \end{equation*}
\begin{equation*}\text{We have } \partial_s T = \alpha N
\text{ for some scalar function } \alpha(t) \in \mathbb{R} \end{equation*}
\begin{equation*}\text{Then } \kappa = \langle \partial_s T, N \rangle \end{equation*}
\begin{equation*}= \langle \alpha N, N \rangle \end{equation*}
\begin{equation*}= \alpha \langle N, N \rangle \end{equation*}
\begin{equation*}= \alpha \|N\|^2 \end{equation*}
\begin{equation*} = \alpha \end{equation*}
i.e. $\alpha = \kappa$
\begin{equation*} \therefore \partial_s T = \alpha N = \kappa N \quad \blacksquare \end{equation*}


2. $ \partial_s N = \kappa T $
\begin{equation*} \text{Since } \|N\| \equiv 1 \end{equation*}
\begin{equation*} \partial_s N \perp N \text{ by the lemma} \end{equation*}
\begin{equation*} \text{We have } \partial_s N = \beta T
\text{ for some scalar function } \beta(t) \in \mathbb{R} \end{equation*}
To determine $\beta$, we differentiate $\langle T, N \rangle = 0$
\begin{equation*} 0 = \partial_s \langle T, N \rangle \end{equation*}
\begin{equation*} = \langle \partial_s T, N \rangle + \langle T, \partial_s N \rangle \end{equation*}
Since $ \partial_s T = \kappa N$ by [Frenet-Serret Theorem](#frenet-serret-theorem) and $\partial_s N = \beta T$
\begin{equation*}
\begin{split}
0 &= \langle \kappa N, N \rangle + \langle T, \beta T \rangle \\
&= \kappa \|N\|^2 + \beta\|T\|^2 \\
&=  \kappa + \beta \quad \text{ as } \|N\|=\|T\|=1 \\
\end{split}
\end{equation*}
\begin{equation*}\therefore \beta = - \kappa \end{equation*}
\begin{equation*} \text{So } \partial_s N = \beta T = - \kappa T \quad \blacksquare \end{equation*}
:::
