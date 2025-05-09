# Evolution Equations

:::{prf:lemma}
:label: commutator

The commutator of arc length and time derivatives is
$$
[\partial_t, \partial_s] = - \kappa^2 \partial_s
$$
:::

- Area
- Length
- Tangent
- Normal
- Curvature


:::{prf:lemma}
:label: evolution-arc-length

The element of arc length $ds$, evolves by:
$$
\partial_t ds = 
$$
:::


:::{prf:proof}

For a single parameter:
Knowing $v = v(u)$
$$ds = v du$$
$$ds = |c'(u,t)| du$$
where $c'$ refers to $\partial_t c$

Differentiate both sides with respect to time:
$$
\partial_t ds = \partial_t (|c'(u)|) du
= \partial_t \partial_u c
= \partial_u \partial_t c
= \partial_u (-\kappa N)
$$


:::
