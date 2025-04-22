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
```

:::{prf:definition}
:label: csf

The Curve Shortening Flow (**CSF**) is the initial value problem
$$
\begin{cases}
\partial_t c &= \vec{\kappa} = -\kappa N \\
c(\cdot, 0) &= c_0 (\cdot)
\end{cases}
$$
:::

Let's take a look at an example.

:::{prf:example}The Shrinking Circle
:label: shrinking-circle

One of the simplest examples of a solution to the CSF is the shrinking circle:
$$
c(u, t) = \sqrt{r_0^2 - 2t} \> (\cos u, \sin u)
$$
where \(r_0\) is the radius at time \(t=0\). See [shrinking circle plot](#shrinking_circle).
:::

```{code-cell} ipython3
:tags: [remove-input]
:label: shrinking_circle

t = np.linspace(0, 2*np.pi, 100)

xval = np.cos(t)
yval = np.sin(t)

fig = px.line(x = xval, y = yval)
fig.update_yaxes(scaleanchor = "x", scaleratio = 1)

# Hide axes
fig.update_layout(xaxis=dict(visible=False), yaxis=dict(visible=False))

# Make plots transparent
fig.update_layout(paper_bgcolor = 'rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')

# Title
fig.update_layout(title_text="Shrinking circle solution of CSF")

fig
```

:::{prf:lemma}
:label: parametric_csf

In a parametrisation, the CSF is the evolution equation
$$
\partial_t (x, y) = d
$$
:::
