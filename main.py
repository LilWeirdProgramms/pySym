import numpy as np
from scipy.integrate import odeint
from bokeh.plotting import figure
from bokeh.layouts import row, column
from bokeh.io import curdoc
from bokeh.models import ColumnDataSource, Slider, Div
from bokeh.models.widgets import Paragraph


"""
To run: bokeh serve --show main.py
"""

rho0 = 1.225
mass = 60
c_w = 1.33
max_area = 10
max_force = {"max": [0, 0]}
time_to_open = 4

args = (mass, c_w, max_area, max_force, time_to_open)


def calc_area(time, args):
    if time < args[4]:
        area = args[2] * time / args[4]
    else:
        area = args[2]
    return area

def height(x, t, *args):
    max_force = args[3]
    height, velocity = x
    grav_const = 9.807 * ((6371 * 1000) / (6371 * 1000 + height)) ** 2
    if height > 0:

        acceleration = args[1] * velocity**2 * 0.5 * \
        rho0 * (288.15 / (288.15 - 0.0065 * height)) ** (1 + ((grav_const * 0.029) / (-8.3314 * 0.0065)))\
        * calc_area(t, args) / args[0] - grav_const

        force = (acceleration + grav_const) * args[0]
        if force > max_force["max"][0]:
            max_force["max"] = [force, t]

        f = [velocity, acceleration]
    else:
        f = [0, 0]
    return f


x0 = [2000, 0]
t = np.linspace(0, 120, 1000)
y = odeint(height, x0, t, args=args)

x1 = [ent/60 for ent in t]
y1 = [entry[0] for entry in y]
y2 = [entry[1] for entry in y]
source1 = ColumnDataSource(data=dict(x=x1, y=y1))
source2 = ColumnDataSource(data=dict(x=x1, y=y2))

TOOLTIP2 = [("Velocity", "@y"), ("Time", "@x"),]
TOOLTIP1 = [("Height", "@y"), ("Time", "@x"),]
fig = figure(plot_height=500, plot_width=800, title="Height over Time:",
              tools="reset,wheel_zoom,hover", tooltips=TOOLTIP1,
              x_range=[0, t[-1]/60], y_range=[0, x0[0] * 1.1])
fig.line("x", "y", source=source1, line_width=3, line_color="red")
fig2 = figure(plot_height=500, plot_width=800, title="Velocity over Time:",
              tools="reset,wheel_zoom,hover", tooltips=TOOLTIP2,
              x_range=[0, t[-1]/60], y_range=[0, min(y2) * 1.1])
fig2.line("x", "y", source=source2, line_width=3)
text_out = Paragraph(text="The Maximum Force which acts on the Fallschirm is %d N, at time = %d s" %(max_force["max"][0],max_force["max"][1]))
drag_slid = Slider(title="Set C_w Fallschirm [-]", value=c_w, start=0.1, end=3, step=0.1)
mass_slid = Slider(title="Set Mass Rocket [kg]", value=mass, start=10, end=1000, step=2)
area_slid = Slider(title="Set Maximum Area Fallschirm [m^2]", value=max_area, start=1, end=100, step=1)
time_slid = Slider(title="Set Simulation Time [s]", value=t[-1], start=2, end=2000, step=2)
height_slid = Slider(title="Set Height in which Fallschirm is opened [m]", value=x0[0], start=100, end=20000, step=100)
vel_slid = Slider(title="Set absolut starting velocity at which Fallschirm is opened [m/s]", value=x0[1], start=0, end=1000, step=5)
open_slid = Slider(title="Set time Fallschirm needs to open [s]", value=time_to_open, start=0.1, end=20, step=0.1)

fig.xaxis.axis_label = "Time in [min]"
fig.yaxis.axis_label = "Height in [m]"
fig2.xaxis.axis_label = "Time in [min]"
fig2.yaxis.axis_label = "Velocity in [m/s]"
fig.title.text_font_size = '16pt'
fig.xaxis.axis_label_text_font_size = "12pt"
fig.yaxis.axis_label_text_font_size = "12pt"
fig2.title.text_font_size = '16pt'
fig2.xaxis.axis_label_text_font_size = "12pt"
fig2.yaxis.axis_label_text_font_size = "12pt"

def update_data(attrname, old, new):
    max_force["max"] = [0, 0]
    args = (mass_slid.value, drag_slid.value, area_slid.value, max_force, open_slid.value)
    t = np.linspace(0, time_slid.value, 1000)
    x0 = [height_slid.value, -vel_slid.value]
    y = odeint(height, x0, t, args=args)
    x1 = [ent / 60 for ent in t]
    y1 = [entry[0] for entry in y]
    y2 = [entry[1] for entry in y]
    text_out.text = "The Maximum Force which acts on the Fallschirm is %d N, at time = %d s" % (max_force["max"][0], max_force["max"][1])
    fig.x_range.start, fig.x_range.end = [0, t[-1]/60]
    fig.y_range.start, fig.y_range.end = [0, x0[0] * 1.1]
    fig2.x_range.start, fig2.x_range.end = [0, t[-1]/60]
    fig2.y_range.start, fig2.y_range.end = [0, min(y2) * 1.1]
    source1.data = dict(x=x1, y=y1)
    source2.data = dict(x=x1, y=y2)


for w in [drag_slid, mass_slid, area_slid, time_slid, height_slid, open_slid, vel_slid]:
    w.on_change('value', update_data)

fig3 = figure(x_range=(0,500), y_range=(0,500), tools="")
fig3.image_url(url=["https://i.ibb.co/rxCwc0P/me.jpg"], x=1, y=1, w=500, h=300, anchor="bottom_left")

fig3.xaxis.visible = False
fig3.yaxis.visible = False
fig3.xgrid.visible = False
fig3.ygrid.visible = False

curdoc().add_root(column(row(fig, fig2), drag_slid, mass_slid, area_slid, time_slid, height_slid, open_slid, vel_slid, text_out, fig3))
curdoc().title = "My Free Fall Calculation"
