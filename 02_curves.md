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

:::{prf: proof}
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

:::{prf: proof}
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