import numpy as np
import plotly.express as px
from plotly.subplots import make_subplots
from plotly import graph_objs as go

# class Planet:
#     def __init__(self, name, radius, eccentricity, semi_minor, semi_major, period, inclination, parent):
#         self.eccentricity = eccentricity
#         self.semi_minor = semi_minor
#         self.semi_major = semi_major
#         self.period = period
#         self.inclination = inclination
#         self.parent = parent
#         self.name = name
#         self.radius = radius
#         self.angle = 0
#         # self.position = Vector2(0, 0)
from typing import NamedTuple
from csv import DictReader


class Planet(NamedTuple):
    name: str
    period: float | int  # in earth days
    semi_major: float | int  # in metres
    eccentricity: float | int  # dimensionless


# AU in metres
AU = 149597870700
# challenge 1:graph orbital period against semi major axis

planets = [
    Planet("mercury", 88, 57.9e6, 0.21),
    Planet("venus", 225, 108.2e6, 0.01),
    Planet("earth", 365.25, 149.6e6, 0.02),
    Planet("mars", 687, 227.9e6, 0.09),
    Planet("jupiter", 4332, 778.6e6, 0.05),
    Planet("saturn", 10759, 1433.5e6, 0.06),
    Planet("uranus", 30685, 2873e6, 0.05),
    Planet("neptune", 60190, 4495.1e6, 0.01),
    Planet("pluto", 90560, 5906.4e6, 0.25),
]

periods = [planet.period / 365 for planet in planets]
semi_majors = [planet.semi_major / AU for planet in planets]
corrected_semi_majors = [pow(planet.semi_major / AU, 1.5) for planet in planets]
eccentricities = [planet.eccentricity for planet in planets]
names = [planet.name for planet in planets]
# pluto orbit angle vs time
from math import cos
import scipy
import numpy as np
from math import pi

pointcnt = 40_000
max_theta = 20
scaler = pointcnt / max_theta
stepsize = 100
for name, period, eccentricity in zip(names, periods, eccentricities):
    x = [x / scaler for x in range(0, pointcnt)]  # from 0 to 20, step size .001
    ys = [pow(1 - cos(theta) * eccentricity, -2) for theta in x]
    I_trap = scipy.integrate.cumulative_trapezoid(ys, x, initial=0)
    I_simp = np.cumsum(
        [
            scipy.integrate.simpson(ys[i : i + stepsize], x[i : i + stepsize])
            for i in range(0, pointcnt - stepsize, stepsize)
        ]
    )
    scaled_I_trap = [
        period * pow(1 - eccentricity * eccentricity, 1.5) * I / (2 * pi)
        for I in I_trap[stepsize::stepsize]
    ]
    scaled_I_simp = [
        period * pow(1 - eccentricity * eccentricity, 1.5) * I / (2 * pi)
        for I in I_simp
    ]
    # swapped x and y axes to be task-accurate
    scaled_xs = [i * period / (2 * pi) for i in x[stepsize::stepsize]]
    import pandas as pd

    data = pd.DataFrame(
        data={
            "Trapeziod method": scaled_I_trap,
            "Simpson's method": scaled_I_simp,
            "zero e approximation": scaled_xs,
        }
    )
    x_degrees = [i * (180 / pi) for i in x]
    fig = px.line(data, y=x[stepsize::stepsize], x=data.columns, title=name)
    fig.update_yaxes(title_text="Angle (Radians)")
    fig.update_xaxes(title_text="Time (Earth Years)")
    fig.show()
    fig.write_html(f"{name}.html")
    # plot I_trap,I_simp, against x
with open("exoplanet.eu_catalog.csv") as exocsv:
    reader = DictReader(exocsv)
    exo_semi_majors = []
    exo_periods = []
    for planet in reader:
        if planet["semi_major_axis"] and planet["orbital_period"]:
            exo_semi_majors.append(planet["semi_major_axis"])
            exo_periods.append(planet["orbital_period"])
    exo_semi_majors = [float(x) for x in exo_semi_majors]
    corrected_exo_semi_majors = [pow(x, 1.5) for x in exo_semi_majors]

fig = make_subplots(rows=1, cols=2, subplot_titles=("Solar System", "Exoplanets"))

exo_scatter = px.scatter(
    x=exo_semi_majors,
    y=exo_periods,
    trendline="ols",
    trendline_options=dict(log_x=True, log_y=True),
    title="Kepler's third law",
)


solar_scatter = px.scatter(
    x=corrected_semi_majors,
    y=periods,
    text=names,
    trendline="ols",
    title="Kepler's third law",

)


exo_traces = []
solar_traces = []
for trace in range(len(exo_scatter["data"])):
    exo_traces.append(exo_scatter["data"][trace])
for trace in range(len(solar_scatter["data"])):
    solar_traces.append(solar_scatter["data"][trace])

for trace in exo_traces:
    fig.append_trace(trace, row=1, col=1)
for trace in solar_traces:
    fig.append_trace(trace, row=1, col=2)

fig.update_xaxes(
    title_text="(semi-major axis in AU) ^1.5 (logarithmic)", type="log", row=1, col=2
)
fig.update_yaxes(
    title_text="Orbit Period in Earth days (logarithmic)", type="log", row=1, col=2
)

solar_scatter.update_xaxes(
    title_text="(semi-major axis in AU) ^1.5 (logarithmic)", type="log"
)
solar_scatter.update_yaxes(
    title_text="Orbit Period in Earth days (logarithmic)", type="log"
)

fig.update_xaxes(
    title_text="(semi-major axis in AU) ^1.5 (logarithmic)", type="log", row=1, col=1
)
fig.update_yaxes(
    title_text="Orbit Period in Earth days (logarithmic)", type="log", row=1, col=1
)
fig.update_traces(textposition="top center", row=1, col=1)
fig.update_traces(textposition="top center", row=1, col=2)

solar_scatter.update_layout(
    updatemenus=[
        dict(
            type="buttons",
            direction="left",
            active=0,
            buttons=list(
                [
                    dict(
                        label="Linear",
                        method="relayout",
                        args=[
                            {
                                "xaxis": {"type": "linear","title_text":"(semi-major axis in AU) ^1.5 (linear)"},
                                "yaxis": {"type": "linear","title_text":"Orbit Period in Earth Days"},
                            }
                        ],
                    ),
                    dict(
                        label="Log",
                        method="relayout",
                        args=[
                            {
                                "xaxis": {"type": "log","title_text":"(semi-major axis in AU) ^1.5 (log)"},
                                "yaxis": {"type": "log","title_text":"Orbit Period in Earth Days"},
                            }
                        ],
                    ),
                ]
            ),
        )
    ]
)
solar_scatter.show()
fig.show()
fig.write_html("plotly_out.html")
print(np.corrcoef(periods, corrected_semi_majors)[0, 1])  # r^2 = 0.999999498!

for planet in planets:
    print(planet.name, (planet.period / 365) / pow((planet.semi_major / AU), 1.5))
