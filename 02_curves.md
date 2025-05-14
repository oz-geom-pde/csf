# Curve Theory

- Parametrisation
- Arc-length
- Length
- Simple, closed curves
- Area
- Curvature
- Frenet-Serret Formulae

# Parametrisation
## Definition of a Parametrised Curve
:::{prf:definition}
:label: parametrised-curve

A parametrised curve is a pair of functions $x(t), y(t)$ that describe a curve in the $\mathbb{R}^2$ plane as:
$$ c(t) = (x(t), y(t)) $$
where $t \in \mathbb{R}$.
:::

## Definition of a Regular Parametrisation
:::{prf:definition}
:label: regular-parametrisation

A curve is called a regular curve if the velocity vector $c'(t)$ is never zero:
$$ c'(t) \neq 0 \quad \text{for all } t $$

i.e. A curve is regular if c(t) is continuous.
:::

# Arc Length
## Arc Length Parameter
:::{prf:definition}
:label: arc-length-parameter

The arc length parameter $L$ measures the distance travelled along a curve from a fixed starting point. 

Given a regular curve $c(t)$, the arc length from $t = a$ to $t$ is:

$$
L(t) = \int_a^t ||c'(u)|| \, du
$$
:::

## Arc Length Parametrisation
:::{prf:definition}
:label: arc-length-parametrisation


:::

## Arc Length of a Segment
:::{prf:definition}
:label: arc-length-segment

The arc length of a segment of a curve from $t = a$ to $t = b$ where $a,b \in \mathbb{R}$ is given by:

$$
L = \int_a^b ||c'(t)|| \, dt
$$
:::

## Total Arc Length
:::{prf:definition}
:label: total-arc-length

The total arc length of a curve $c(t)$ over the full interval of $c(t): t \in [a, b]$ is:

$$ L = \int_a^b ||c'(t)|| \, dt $$

This gives the total distance travelled along the curve from start to end.
:::

# Simple, Closed Curves
## Definition of Simple, Closed Curves
:::{prf:definition}
:label: simple-closed-curve

A curve $c: [a, b] \in \mathbb{R}^2$ is called a simple closed curve if:

- It is smooth on $[a, b]$: $c$ is continuously differentiable on $[a, b]$
- It starts and ends at the same point: $c(a) = c(b)$
- It is injective on $(a, b)$: $c(t_1) \ne c(t_2)$ for all $a < t_1 < t_2 < b$
- The endpoints match smoothly: $c'(a) = c'(b)$

This means the curve traces out a loop without crossing itself, and joins up smoothly at the endpoints. 
:::

# Curvature
## Definition of Curvature
:::{prf:definition}
:label: curvature-definition

Let $c(t) = (x(t), y(t))$ be a regular smooth curve parametrised by arc-length $s$. The **curvature** $\kappa(s)$ at a point on the curve is defined as:
$$ \kappa(s) = \langle \partial_s T, N \rangle $$

Here:
- $T = \partial_s c$ is the unit tangent vector.
- $N = R_{\pi/2}$ is the unit normal vector.

This gives $\kappa$ a signed value, depending on whether the curve bends to the left (positive) or right (negative).
:::

:::{prf:remark}
:label: curvature-interpretation

Curvature measures the rate at which the tangent vector turns as a particle moves along the curve. It captures how intense the "bend" of the curve is at each point:
- If $\kappa(s) = 0$, the curve is straight.
- As $||\kappa(s)||$ increases, a sharper curve is indicated at that particular point.
:::

## Dependence on Unit Normal

:::{prf:remark}
:label: unit-normal-dependence

The curvature vector $\partial_s T$ points in the direction of the unit normal vector $N$, i.e.:

$$
\partial_s T = \kappa N
$$

However, since $N$ is a general case, and broadly defined, both $N$ and $-N$ are orthogonal to $T$, this introduces a **sign ambiguity** in the curvature.

To resolve this, we've chosen the orientation of the plane so that counterclockwise is positive. From here, we defined $N = R_{\pi/2}$ so that the pair $(T, N)$ is **positively oriented** (i.e. turning from $T$ to $N$ is counter-clockwise). Therefore $\kappa$ is well-defined with a consistent sign/direction.
:::



# Frenet-Serret Formulae
## Product Rule
:::{prf:lemma}
:label: Product rule

$$\text{If } X, Y: (a,b) \rightarrow \mathbb{R}^2 $$
Where $X(t) = (x_1(t), x_2(t))$ and $Y(t) = (y_1(t), y_2(t))$
$$
\text{Then } \frac{d}{dt} \langle X(t), Y(t) \rangle = 
\langle \frac{d}{dt} X(t), Y(t) \rangle + \langle X(t), \frac{d}{dt} Y(t) \rangle
$$

$$
\text{Using dot product :}
\quad \frac{d}{dt}(X \cdot Y) = \frac{d}{dt}X \cdot Y + X \cdot \frac{d}{dt}Y
$$

$$
(X \cdot Y)' = X' \cdot Y + X \cdot Y'
$$
:::

:::{prf:proof}
:label: Product rule
$$ (X \cdot Y)' = (x_1 y_1 + x_2 y_2)' $$
$$ = (x_1 y_1)' + (x_2 y_2)' $$
$$ = x_1' y_1 + x_1 y_1' + x_2' y_2 + x_2 y_2' $$
$$ = (x_1', x_2') \cdot (y_1, y_2) + (x_1, x_2) \cdot (y_1', y_2') $$
$$ \therefore = X' \cdot Y + X \cdot Y' \quad \blacksquare $$
:::

## Orthogonality of the derivative of a constant length vector field

:::{prf:lemma}
:label: 

$$ \text{If } ||V|| \equiv \alpha $$
$$ \text{then } \partial_S V \perp V $$
$$ \text{i.e. } \partial_S V \cdot V = 0 $$
$$ \therefore \text{Similarly } \partial_S T \perp T = 0 $$
$$ \partial_S N \perp N = 0 $$
:::

:::{prf:proof}
:label:

Since $||V(u)|| \equiv 1$ (i.e. for every u, $||{V(u)}|| = 1$)
$$ = \partial_S (1) $$
$$ = \partial_S ||V||^2 $$
$$ = (T \cdot T)' \quad \text{ here, } T' \text{ refers to } \partial_u T $$
$$ = (T' \cdot T + T \cdot T') \text{  by product rule} $$
$$ 0 = 2 (T' \cdot T) $$ 
Then we must have $T' \cdot T = 0$
$$ \therefore \partial_S T \cdot T = \frac{1}{v} T' \cdot T = 0 \quad \blacksquare $$
:::

## Frenet-Seret Equations
:::{prf:theorem}
:label: frenet-serret-theorem

$$ \partial_s T = \kappa N $$
$$ \partial_s N = - \kappa T $$

Or, in matrix form:
$$ \partial_s
\begin{pmatrix}
T \\ N
\end{pmatrix}
= \begin{pmatrix}
0 & \kappa \\ -\kappa & 0
\end{pmatrix}
\begin{pmatrix}
T \\ N
\end{pmatrix}
$$
:::


:::{prf:proof}
:label: frenet-serret-proof
1. $ \partial_s T = \kappa N $
$$ \text{By definition: } \quad \kappa = \langle \partial_s T, N \rangle $$
$$ \text{Since } \partial_s T \perp T $$
$$ \text{We have } \partial_s T = \alpha N
\text{ for some scalar function } \alpha(t) \in \mathbb{R} $$
$$ \text{Then } \kappa = \langle \partial_s T, N \rangle $$
$$ = \langle \alpha N, N \rangle $$
$$ = \alpha \langle N, N \rangle $$
$$ = \alpha ||N||^2 $$
$$ = \alpha $$
i.e. $\alpha = \kappa$
$$ \therefore \partial_s T = \alpha N = \kappa N \quad \blacksquare $$


2. $ \partial_s N = \kappa T $
$$ \text{Since } ||N|| \equiv 1 $$
$$ \partial_s N \perp N \text{ by the lemma} $$
$$ \text{We have } \partial_s N = \beta T
\text{ for some scalar function } \beta(t) \in \mathbb{R} $$
To determine $\beta$, we differentiate $\langle T, N \rangle = 0$
$$ 0 = \partial_s \langle T, N \rangle $$
$$ = \langle \partial_s T, N \rangle + \langle T, \partial_s N \rangle $$
Since $ \partial_s T = \kappa N$ by (1) and $\partial_s N = /beta T$ by (2)
$$ 0 = \langle \kappa N, N \rangle + \langle T, \beta T \rangle $$
$$ = \kappa ||N||^2 + \beta||T||^2 $$
$$ = \kappa + \beta \quad \text{ as } ||N||=||T||=1 $$
$$ \therefore \beta = - \kappa $$
$$ \text{So } \partial_s N = \beta T = - \kappa T \quad \blacksquare $$
:::