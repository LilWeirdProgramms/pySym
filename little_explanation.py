from scipy.integrate import odeint
from numpy import linspace
from math import exp

"""
In this File the Differential Equation:
    y'' + y' = 2y       (Solution y = e^t)
is solved.

To solve it, it is transformed to 'Zustandsraumdarstellung':
y = x1
x1' = x2 = y'
x2' = 2x1 - x2
"""


def derivation(y, t):

    """
    needs to be of form: dy/dt = f(y,t)
    In:
        - Initial Conditions for Parameters to be calculated.
        - array of times for which to solve
    Returns:
        Derivates of Parameters
    """

    x1, x2 = y
    dx1_dt = x2
    dx2_dt = 2*x1 - x2

    return [dx1_dt, dx2_dt]


eval_at = 100
t_vec = linspace(0, eval_at, num=100)    # System is solved for 0 <= t <= eval_at with 100 steps
init_values = [1, 1]    # Initial Values: y(0) = 1, y'(0) = 1
res = odeint(derivation, init_values, t_vec)    # Integrates with RK45; returns the values of y for each t

print("Numeric = " + str([exponetial[0] for exponetial in res][-1]), "Anal = " + str(exp(eval_at)))


"""
Plotting:

import bokeh.plotting as plt
from bokeh.models import ColumnDataSource
plt.output_file("this.html")

t_vec2 = linspace(0,10,num = 5)
x0 = [1, 1]
res2 = odeint(deriv, x0, t_vec2)

plot_res = [x[0] for x in res]
plot_res2 = [x[0] for x in res2]
fig = plt.figure()
dici = {"xs":[t_vec, t_vec2], "ys":[plot_res, plot_res2], "color": ["red", "blue"], "line_width":[2,2],
"legend":["100", "not 100"]}
sourcei = ColumnDataSource(dici)
fig.multi_line(xs="xs", ys="ys", source=sourcei, color="color", line_width = "line_width", legend= "legend")
plt.show(fig)
"""