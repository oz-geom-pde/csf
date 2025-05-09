# Curve Theory

- Parametrisation
- Arc-length
- Length
- Area
- Curvature

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