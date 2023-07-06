import numpy as np


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


class Planet(NamedTuple):
    name: str
    period: float | int  # in earth days
    semi_major: float | int  # in metres


# AU in metres
AU = 149597870700
# challenge 1:graph orbital period against semi major axis

planets = [Planet("mercury", 88, 57.9e6), Planet("venus", 225, 108.2e6),
           Planet("earth", 365.25, 149.6e6), Planet("mars", 687, 227.9e6),
           Planet("jupiter", 4332, 778.6e6), Planet("saturn", 10759, 1433.5e6),
           Planet("uranus", 30685, 2873e6), Planet("neptune", 60190, 4495.1e6),
           Planet("pluto", 90560, 5906.4e6)]

periods = [planet.period/365 for planet in planets]
corrected_semi_majors = [pow(planet.semi_major/AU, 1.5)
                         for planet in planets]
print(np.corrcoef(periods, corrected_semi_majors)
      [0, 1])  # r^2 = 0.999999498!

for planet in planets:
    print(planet.name, (planet.period / 365) /
          pow((planet.semi_major / AU), 1.5))
