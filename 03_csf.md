---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Curve Shortening Flow

```{code-cell} ipython3
:tags: [remove-input]

import numpy as np
import sympy as sp

import pandas as pd

import plotly.express as px
import plotly.graph_objects as go

import plotly.io as pio
pio.templates.default  = "plotly_dark"

import imageio.v3 as iio

import lib.flow_plot as fp
```

## Definition of the curve shortening flow

:::{tip} Definition
:icon: false
:label: csf
Let $\gamma_0 : I \rightarrow \mathbb{R}^2$ be a regular parametrized curve. **The curve shortening flow (CSF)** is a family of regular parametrized curves $\gamma_t : I \rightarrow \mathbb{R}^2$ for $t \in \left[0, T \right >$ such that the function $\gamma : \left[0, T \right > \times I \rightarrow \mathbb{R}^2$ defined by $\gamma(t,u) = \gamma_t(u)$ is a smooth solution to the initial value problem
\begin{equation*}
\begin{cases}
\partial_t \gamma = -\vec{\kappa} = - \kappa N \\
\gamma(0, \cdot) = \gamma_0(\cdot)
\end{cases}
\end{equation*}
where $\kappa(t,\cdot), N(t, \cdot)$ are the curvature and the unit normal of $\gamma_t$ respectively.
:::

:::{prf:remark}
:numbered: false
Both the curvature $\kappa$ and the unit normal $N$ depend on the orientiation of $\gamma$, in the sense that changing the orientation of the parametrization changes the signs of both $\kappa$ and $N$. Therefore, the curvature vector $\vec{\kappa}$ does not depend on the choice of orientation of $\gamma.$
:::

## Invariance properties of the curve shortening flow

Before we consider explicit examples of the curve shortening flow, we prove two important facts about the curve shortening flow:
- the curve shortening flow is invariant under isometries of $\mathbb{R}^2$
- the curve shortening flow is invariant under tangential perturbations

:::{danger} Proposition
:icon: false
Let $\gamma_0 : I \rightarrow \mathbb{R}^2$ be a regular parametrized curve, let $\gamma$ be its curve shortening flow and let $A : \mathbb{R}^2 \rightarrow \mathbb{R}^2$ be an isometry of $\mathbb{R}^2$. <br>
Then $A \circ \gamma$ is the curve shortening flow of $A \circ \gamma_0$.
:::

:::{prf:proof}
:numbered: false
Since $A$ is an isometry of $\mathbb{R}^2$, it is in particular a linear map $A \in L(\mathbb{R}^2, \mathbb{R}^2).$ Denote by $A$ its matrix representation in the canoninal basis $\{e_1, e_2\}$ of $\mathbb{R}^2$ so that 
\begin{equation*}
A(x,y) = A \cdot \begin{bmatrix} x \\ y \end{bmatrix}.
\end{equation*}
By direct computation, we have:
\begin{equation*}
\begin{split}
\frac{\partial}{\partial t} \left( A \circ \gamma \right) &= \frac{\partial}{\partial t} \left(A \cdot \gamma \right) \\
&= A \cdot \frac{\partial \gamma}{\partial t} \\
&= A \cdot \left( - \kappa N \right) \\
&= -\kappa \cdot A \cdot N \\
&= - \kappa \cdot (A \circ N) \\
&= - \kappa_{A \circ \gamma} \cdot N_{A \circ \gamma}.
\end{split}
\end{equation*}
Since $(A \circ \gamma)(0, \cdot) = A(\gamma(0,\cdot)) = A(\gamma_0(\cdot)) = (A \circ \gamma_0 )(\cdot)$, we see that $A \circ \gamma$ is the curve shortening flow of $A \circ \gamma_0.$
:::

:::{danger} Proposition
:icon: false
:label: prop_reparam
Let $\gamma_0 : I \rightarrow \mathbb{R}^2$ be a regular parametrized curve, let $\gamma$ be its curve shortening flow and let $\varphi: J \rightarrow I$ be an orientation-preserving diffeomorphism. <br>
Then $\tilde{\gamma}(t,v) \coloneqq \gamma(t,\varphi(v))$ is the curve shortening flow of the reparametrization $\gamma_0 \circ \varphi.$
:::

:::{prf:proof}
:numbered:false
By direct computation, we have:
\begin{equation*}
\frac{\partial \tilde{\gamma}}{\partial t}(t,v) &= \frac{\partial}{\partial t} (\gamma(t,\varphi(v))) \\
&=  D \gamma(t,\varphi(v)) \cdot \begin{bmatrix}1 \\ 0 \end{bmatrix} \\
&= \frac{\partial \gamma}{\partial t}(t,u) \\
&= - \kappa(t,u) N(t,u) \\
&= - \kappa(t, \varphi(v)) N(t, \varphi(v)).
\end{equation*}
Since $\tilde{\gamma}(0,v) = \gamma(0,\varphi(v)) = \gamma_0(\varphi(v)) = (\gamma_0 \circ \varphi)(v)$, we see that $\tilde{\gamma}$ is the curve shortening flow of $\gamma_0 \circ \varphi.$
:::

:::{prf:remark}
:numbered: false
:label: rem_tangper
From the [previous lemma](#prop_reparam) we know that the curve shortening flow is invariant under reparametrizations of the initial curve, meaning that the way in which the initial curve is parametrized does not play a role in the existence of the solution to the flow. Furthermore, there is even a stronger invariance of the curve shortetning flow - **the invariance under tangential perturbations.** Namely, if the family of parametrized curves $\gamma : \left[ 0, T \right> \times I \rightarrow \mathbb{R}^2$ solves the initial value problem
\begin{equation*}
\begin{cases}
\displaystyle \frac{\partial \gamma}{\partial t}(t,u) = -\kappa(t,u)N(t,u) + \alpha(t,u)T(t,u) \\
\gamma(0, \cdot) = \gamma_0(\cdot)
\end{cases}
\end{equation*}
then there exists a family of time-dependent diffeomorphisms $\varphi : \left[0,T\right> \times J \rightarrow I$ such that $\varphi(0, \cdot) = \text{id}(\cdot)$ and 
\begin{equation*}
\frac{\partial}{\partial t} \gamma(t, \varphi(t,u)) = -\kappa(t,\varphi(t,u))N(t,\varphi(t,u)),
\end{equation*}
meaning that the family $\gamma(t, \varphi(t,u))$ is the curve shortening flow of
\begin{equation*}
\gamma(0,\varphi(0,\cdot)) = \gamma_0(\varphi(0,\cdot)) = \gamma_0(\cdot).
\end{equation*}
In other words, it is not necessary to require that the solution of the curve shortening flow has a vanishing tangential component, since such a tangential component arises only due to a particular choice of parametrization of the initial curve.
:::

:::{note} Corollary
:icon: false
Let $\gamma_0 : I \rightarrow \mathbb{R}^2$ be a regular parametrized curve and let $\gamma :  \left[ 0, T \right > \times I \rightarrow \mathbb{R}^2$ be a smooth family of parametrized curves such that $\gamma(0, \cdot) = \gamma_0$ and that satisfies
\begin{equation*}
\langle \frac{\partial \gamma}{\partial t}, N \rangle = -\kappa.
\end{equation*}
Then, such a family can locally be reparametrized to a curve shortening flow. If $I$ is compact, such a reparametrization can be achieved globally.
:::

Because of these invariance properties, when we talk about the curve shortening flow we can consider curves in question simply as subsets of $\mathbb{R}^2$, without referencing their parametrizations.

## Examples
Let's take a look at some examples.

:::{caution} Example (The Shrinking Circle)
:icon:false
:label:exa_circle
Let $r_0 > 0$ be a positive number and let $\gamma_0 : \mathbb{R} \rightarrow \mathbb{R}^2$ be a parametrization of the circle with centre $(0,0)$ and radius $r_0$ given by
\begin{equation*}
\gamma_0(u) = r_0(\cos u, \sin u).
\end{equation*}
Keeping in mind the geometric aspect of the curve shortening flow, a reasonable guess for the solution of the curve shortening flow of $\gamma_0$ is the family 
\begin{equation*}
\gamma(t,u) = r(t) \cdot \gamma_0(u).
\end{equation*}
By direct computation if follows that
\begin{equation*}
\kappa(t,u) = \frac{1}{r_0 \cdot r(t)}, \quad N(t,u) = \gamma_0(u)
\end{equation*}
so the curve shortening flow boils down to the following initial value problem
\begin{equation*}
\begin{cases}
\displaystyle r' = - \frac{1}{r} \\
r(0) = r_0
\end{cases}
\end{equation*}
It is straight-forward to check that the solution is given by $r(t) = \sqrt{r_0^2 - 2t}$ so the curve shortening flow of $\gamma_0$ is given by
\begin{equation*}
\gamma(t,u) = \sqrt{r_0^2-2t} \cdot (\cos u, \sin u).
\end{equation*}
Therefore, under the curve shortening flow, the circle shrinks to a point in finite time.
:::

:::{prf:remark}
:numbered: false
From the [example above](#exa_circle) it follows that the maximal time of existence of the curve shortening flow of the circle with initial radius $r_0$ is equal to
\begin{equation*}
T_{\text{max}} = \frac{r_0^2}{2}.
\end{equation*}
:::

:::{figure} ./video/circle.mp4
The curve shortening flow of the circle with initial radius $r_0 = 2$. Under the flow, the circle shrinks to a point in finite time that equals $T_{\text{max}} = 2.$
:::

:::{caution} Example
:icon: false
Consider the family of parametrized curves given by 
\begin{equation*}
\gamma(t,u) = \sqrt{- 2t} \cdot (\cos(u + e^t), \sin(u+e^t)).
\end{equation*}
By direct computation it follows that
\begin{equation*}
\begin{split}
T(t,u) &= (-\sin(u+e^t), \cos(u+e^t)), \quad N(t,u) = - (\cos(u+e^t), \sin(u+e^t)) \\
\kappa(t,u) &= \frac{1}{\sqrt{-2t}}
\end{split}
\end{equation*}
and so
\begin{equation*}
\begin{split}
\frac{\partial \gamma}{\partial t}(t,u) &= - \frac{1}{\sqrt{-2t}}(\cos(u+e^t), \sin(u+e^t)) + e^t \sqrt{-2t} (-\sin(u+e^t), \cos(u+e^t)) \\
&= - \kappa N + e^t \sqrt{-2t} \cdot T.
\end{split}
\end{equation*}
Because the curve shortening flow is [invariant under tangential perturbations](#rem_tangper), it follows that the family $\gamma(t,u)$ is the curve shortening flow of some curve. Since $\gamma(0,\cdot)$ parametrizes a unit circle $\mathbb{S}^1,$ the family $\gamma(t,u)$ is precisely the curve shortening flow of the unit circle. This does not contradict the computations made in the [previous example](#exa_circle) because we can achieve one solution from the other by reparametrization.
:::


:::{caution} Example (The Grim Reaper)
:icon: false
Let $f \colon \langle -\pi/2, \pi / 2 \rangle \rightarrow \mathbb{R}$ be defined by $f(x) \coloneqq - \ln \left(\cos x \right)$. The graph of the function $f$ is a curve that can be parametrized by
\begin{equation*}
\gamma_0(u) = (u, f(u)) = (u, -\ln\left(\cos u \right)), \quad u \in \langle -\pi/2, \pi / 2 \rangle.
\end{equation*}
We call $\gamma_0$ the **Grim Reaper.** <br>
A straight-forward computation shows that the family of curves $\gamma = \gamma(t,u)$ given by 
\begin{equation*}
\gamma(t,u) = (u, f(u) + t) = (u, - \ln\left( \cos u \right) + t)
\end{equation*}
is the curve shortening flow of $\gamma_0$. Note that $\displaystyle \frac{\partial \gamma}{\partial t} = (0,1)$ so under the curve shortening flow $\gamma_0$ moves by vertical translations in the direction of $y$-axis.
:::

<!--```{code-cell} ipython3
:tags: [remove-input]
:label: shrinking_circle

uvals = np.linspace(0, 2*np.pi, 200)
r0 = 3
tvals = np.linspace(0, 0.5*(r0**2), 20)

u = fp.u
t = fp.t

r = (r0**2 - 2*t)**0.5
xfunc = r * sp.cos(u)
yfunc = r * sp.sin(u)

flow = fp.ExplicitFlow(uvals = uvals, tvals = tvals, xfunc = xfunc, yfunc = yfunc)

# Hide axes
flow.fig.update_layout(xaxis=dict(visible=False), yaxis=dict(visible=False))

# Make plots transparent
flow.fig.update_layout(paper_bgcolor = 'rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')

# Title
flow.fig.update_layout(title_text="Shrinking circle solution of CSF")
flow.fig
```-->

## Parabolicity
In this section, our goal is to understand what kind of an equation is the curve shortening flow.

:::{tip} Definition
:icon: false
Let $U \subseteq \mathbb{R}^n$ be an open set and let $U_T = \left < 0, T \right] \times U$. A second-order equation of the form
\begin{equation*}
u_t + Lu = f,
\end{equation*}
where for each time $L$ is a second-order partial differential operator given by 
\begin{equation*}
Lu = - \sum_{i,j = 1}^{n}a^{ij}(t,x)u_{x_i, x_j} + \sum_{i = 1}^n b^i(t,x)u_{x_i} + c(t,x)u,
\end{equation*}
is **parabolic** if there exists a constant $\theta > 0$ such that
\begin{equation*}
\sum_{i,j = 1}^n a^{ij}(t,x)\xi_i \xi_j \geq \theta \lvert \xi \rvert^2,
\end{equation*}
for all $(t,x) \in U_T, \xi \in \mathbb{R}^n$.
:::
In particular, the matrix of a parabolic system is positive definite.

:::{prf:remark}
:numbered: false
- An obvious example of a parabolic equation is the *heat equation*.
- General second-order parabolic equations describe diffusion processes. In what follows, we will see that the curve shortening flow describes diffusion of curvature.
:::

Now, let $\gamma_0 : I \rightarrow \mathbb{R}^2$ be a regular parametrized curve and let $\gamma = \gamma(t,s)$ be its curve shortening flow, where we can without loss of generality assume that $s$ is the arc-length parameter, meaning that $\gamma(t, \cdot)$ is parametrized by arc-length $s = s(t)$. From the Frenet-Serret equations we know that $\gamma_{ss} = -\kappa N,$ so the curve shortening flow can be written as
\begin{equation*}
\gamma_t = \gamma_{ss}.
\end{equation*}
Written in this form, the curve shortening flow looks like a linear heat equation, however that is hardly the case. <br>
First of all, the arc-length $s = s(t)$ depends nonlinearly on time $t$, so $\gamma_{ss}$ does not depend linearly on $\gamma$. Therefore, the curve shortening flow is not a linear equation. <br>
Secondly, the curve shortening flow is not strictly parabolic. Indeed, if $\gamma = \gamma(t,s)$ is the curve shortening flow, then we have
\begin{equation*}
\begin{split}
\frac{\partial \gamma}{\partial t} &= \frac{\partial^2 \gamma}{\partial s^2} \\
&= \frac{1}{\lvert \gamma' \rvert} \frac{\partial}{\partial u} \left( \frac{1}{\lvert \gamma' \rvert} \frac{\partial \gamma}{\partial u} \right) \\
&= \frac{1}{\lvert \gamma' \rvert} \left( \frac{\partial^2 \gamma}{\partial u^2} - \langle \frac{\gamma'}{\lvert \gamma' \rvert}, \frac{\partial^2 \gamma}{\partial u^2} \rangle \frac{\gamma'}{\lvert \gamma' \rvert} \right)
\end{split}
\end{equation*}
In components $\gamma = (x,y)$, the equation above can be written as 
\begin{equation*}
\begin{bmatrix} x_t \\ y_t \end{bmatrix} = \frac{1}{\lvert \gamma' \rvert^4} \begin{bmatrix} y_u^2 & -x_u y_u \\ -x_u y_u & x_u^2 \end{bmatrix} \cdot \begin{bmatrix} x_{uu} \\ y_{uu} \end{bmatrix}.
\end{equation*}
The determinant of the matrix of this system is equal to $0$ and so it is only positive semidefinite, and not positive definite. Therefore, the curve shortening flow is not strictly parabolic equation, but rather degenarate parabolic equation. 

This degeneracy is due to the fact that the curve shortening flow is [invariant under tangential perturbations](#rem_tangper), and it has important consequences regarding the question on existence and uniqueness of solution to the curve shortening flow given an initial curve $\gamma_0$. Because the curve shortening flow is not strictly parabolic, we cannot use the standard PDE theory cannot be applied to prove existence and uniqueness of its solutions. In order to prove such a result, we need to overcome the degeneracy of the equation of the flow. There are various way in which this can be achieved, one of which is by expressing the flow as an graph over $\gamma_0$, i.e. using the ansatz 
\begin{equation*}
\gamma(t,u) = \gamma_0(u) + f(t,u)N(u),
\end{equation*}
where $f$ is some smooth function. Using this ansatz, the curve shortening flow boils down to a quasilinear strictly parabolic equation on $f$, to which the standard PDE theory can be applied. It is worth noting that using this metod, the solution $\gamma$ defined by the ansatz above only solves the curve shortening flow up to tangential motion, but any tangential component can be removed by a reparametrization.

:::{danger} Theorem (Short-time existence of the curve shortening flow)
:icon: false
Let $\gamma_0 : I \rightarrow \mathbb{R}^2$ be a regular parametrized curve. Then there exists $T > 0$ and a smooth function $\gamma : \left[ 0, T \right> \times I \rightarrow \mathbb{R}^2$ which solves the curve shortening flow
\begin{equation*}
\begin{cases}
\displaystyle \frac{\partial \gamma}{\partial t} = - \kappa N \\
\gamma(0, \cdot) = \gamma_0(\cdot)
\end{cases}
\end{equation*}
:::