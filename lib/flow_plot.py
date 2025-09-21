import numpy as np
import sympy as sp

import pandas as pd

import plotly.express as px
import plotly.graph_objects as go

import plotly.io as pio
pio.templates.default  = "plotly_dark"

u = sp.symbols("u")
t = sp.symbols("t")

class ExplicitFlow:
    def __init__(self, uvals, tvals, xfunc, yfunc):
        U, T = np.meshgrid(uvals, tvals)
        self.U = U.flatten()
        self.T = T.flatten()

        xlambda = sp.lambdify((u, t), xfunc)
        self.X = xlambda(u = self.U, t = self.T)

        ylambda = sp.lambdify((u, t), yfunc)
        self.Y = ylambda(u = self.U, t = self.T)

        self.plot = FlowPlot(self.X, self.Y, self.T)
        self.fig = self.plot.fig

class FlowPlot:
    def __init__(self, x, y, t):
        self.df = pd.DataFrame({"x": x, "y": y, "t": t})

        df0 = self.df.query('t==0.0')
        xmin = df0.x.min()
        xmax = df0.x.max()
        ymin = df0.y.min()
        ymax = df0.y.max()

        pad = 0.1 * (ymax-ymin)

        self.fig = px.line(self.df, x="x", y="y", animation_frame="t", range_x=[xmin-pad,xmax+pad], range_y=[ymin-pad, ymax+pad])

        self.fig.update_yaxes(scaleanchor = "x", scaleratio = 1)
