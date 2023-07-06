#include "raylib.h"
#include "raymath.h"
#include <stdio.h>
#define abs(x) ((x) < 0 ? -(x) : (x))
#define max(x, y) ((x) > (y) ? (x) : (y))
#define min(x, y) ((x) < (y) ? (x) : (y))
// I build with cc nbody_sim.c -O3 -Wall -Wextra -pedantic -ffast-math -funsafe-math-optimizations -l:libraylib.a -lm -pthread -o nbody_sim

// http://www.eclecticon.info/index_htm_files/BPhO%20CompPhys%2006%20Planets.pdf
// https://www.dropbox.com/s/dot9l0igz7x1ija/Comp%20Challenge%20%202023%20Solar%20System%20orbits%20-%20presentation.pdf?dl=0

// https://ssd.jpl.nasa.gov/ has all the solar system data!
// https://exoplanetarchive.ipac.caltech.edu/ has some exoplanets
// challenges:
// challenge 1: Replicate kepler's 3rd law correlation
// challenge 2: plot elliptical orbits of planets
// assume sun is stationary origin, orbits are ellipses in same plane
// challenge 3: create a 2d animation of the planets orbiting the sun
// plot inner 5 planets seperately from outer planets
// inner planets have 1 earth orbit per s, outer planets have 1 jupiter orbit per s
// challenge 4: create a 3d animation of the planets orbiting the sun
// Pluto is basically the only one which will be noticeably inclined, but add data for all!
// challenge 5: calculate orbit angle vs time for an eccentric orbit e.g. Pluto
// and compare to a circular one with same period
// will have to do numerical integration with cumsum, simpson's rule
// challenge 6: solar system spirograph
// Choose a pair of planets and determine their orbits vs time. At time intervals of
// Dt, draw a line between the planets and plot this line. Keep going for N
// orbits of the outermost planet.
// challenge 7: like Ptolemy, use orbital models to plot the orbits of bodies relative to
// another orbiting body, e.g. earth
// extensions!
// make n body sims of binary stars!
// do more funky exoplanet stuff
// data can be found at
// visualise lagrange points
// relativistic corrections
// {\displaystyle \sigma ={\frac {24\pi ^{3}L^{2}}{T^{2}c^{2}\left(1-e^{2}\right)}}\ ,}
// sigma = (24 pi^3 a^2)/(T^2c^2(1-e^2)) where sigma is the perihelion shift
// in radians per revolution (only really relevant for Mercury and Venus)

// https://www.youtube.com/watch?v=Qw6uZQYwG7A
// write up a report/paper in latex
// can do wasm stuff to make a website!
// visualise transfer orbits (hohmann, const thrust (nuclear, solar sails, ion engines)),
// gravity assists, etc.
// variables
// e (normally \epsilon) = eccentricity, a = semi-major axis,
// b = semi-minor axis,T = orbital period, beta = orbital inclination in degrees

// Equations:
// F = GmM/r^2
// F = ma
// T^2 = ((4*pi^2)/G(M+m)) * a^3
// e = sqrt(1-(b^2/a^2))
// polar equation of ellipse where
// r = a(1-e^2)/(1+e*cos(theta))

// kepler's 3 laws:
// 1. the orbit of each planet is an ellipse with the sun at one focus
// 2. a line joining a planet and the sun sweeps out equal areas during equal
// intervals of time
// dA/dt = 1/2sqrt(G(m+M)(1-e^2)/a) (this is constant)
// 3. the square of the orbital period of a planet is proportional to the cube
// of the semi-major axis of its orbit

// notes:
// because (M_sun + m_planet) ≈ M_sun, we can ignore mass of planet, finding the following equation:
// (orbital period/Earth Year) ≈ (a/AU)^3/2
// (This can be graphed for both our solar system and exoplanet.eu data)
// line of best fit can be found and grad diff from linear should be v small
// (remember, best fit should go through 0,0)

// there are various ways to detect exoplanets (see here: https://en.wikipedia.org/wiki/Methods_of_detecting_exoplanets)

// semi major axis can be found by max dist from barycentre
// period can be found by time to start pos within epsilon e.g. (in AU) 1e-6
typedef struct Body {
    // probably should use SI units I guess
    struct Vector3 pos;
    struct Vector3 vel;
    struct Vector3 acc;
    double mass;
    double radius;
    const char *name;
    Color color; // maybe later add texture, atmosphere, etc.
} Body;

typedef Body Planet;

typedef struct myPlanet {
    struct {
        struct Vector3 pos;
        struct Vector3 vel;
        struct Vector3 acc;
        double mass;
        double radius;
        const char *name;
        Color color;
    };
    double semiMajorAxis;
    Vector3 orbitTrail[10000]; // Is this enough?
} myPlanet;

#define BODYCNT 6
int main() {
    InitWindow(800, 450, "raylib planets - basic demo?");
    Vector3 centre = {GetScreenWidth() / 2.0f, GetScreenHeight() / 2.0f, 0};
    SetTargetFPS(60);

    static const float G = 6.67408e-11 * scale * scale;
    static const float scale = 5e-3;
    static const float AU = 1.496e11 * scale;        // delta of 1e4
    static const float DT = 3600.;                   // 1 hour
    static const float TIMELIM = 3600. * 24. * 365.; // 100 year
    float maxXY = 0;
    Body sun = {
        .name = "Sun",
        .pos = {0, 0, 0},
        .vel = {0, 0, 0},
        .acc = {0, 0, 0},
        .mass = 1.989e30,
        .radius = 6.957e8 * 5e-2, // sun is too big!!
        .color = RED,
    };

    Planet earth = {
        .name = "Earth",
        .pos = {AU, 0, 0},
        .vel = {0, 29784, 0},
        .acc = {0, 0, 0},
        .mass = 5.972e24,
        .radius = 6.371e6,
        .color = BLUE,
    };
    Body moon = {
        .name = "Moon",
        .pos = {AU + 3.844e8 * scale, 0, 0},
        .vel = {0, 29784 + 1022, 0},
        .acc = {0, 0, 0},
        .mass = 7.34767309e22,
        .radius = 1.7374e6,
        .color = GRAY,
    };

    Planet jupiter = {
        .name = "Jupiter",
        .pos = {5.204267 * AU, 0, 0},
        .vel = {0, 13070, 0},
        .acc = {0, 0, 0},
        .mass = 1.898e27,
        .radius = 6.9911e7,
        .color = ORANGE,
    };
    Planet saturn = {
        .name = "Saturn",
        .pos = {9.582017 * AU, 0, 0},
        .vel = {0, 9690, 0},
        .acc = {0, 0, 0},
        .mass = 5.683e26,
        .radius = 5.8232e7,
        .color = BROWN,
    };
    Planet mercury = {
        .name = "Mercury",
        .pos = {0.38709893 * AU, 0, 0},
        .vel = {0, 47360, 0},
        .acc = {0, 0, 0},
        .mass = 5.972e24,
        .radius = 6.371e6,
        .color = LIGHTGRAY,
    };

    Body bodies[BODYCNT] = {sun, mercury, earth, moon, jupiter, saturn};
    Planet *planets = bodies; // can exclude sun as static iguess
    for (int i = 0; i < BODYCNT; i++) {
        bodies[i].vel = Vector3Scale(bodies[i].vel, scale * scale);
        bodies[i].mass = bodies[i].mass * scale * scale * scale;
        bodies[i].radius = bodies[i].radius * scale;
        printf("%s: mass: %f, radius: %f\n", bodies[i].name, bodies[i].mass, bodies[i].radius);
    }
    float t = 0;
    size_t iter = 0;
    const int subLim = 1000; // sub iterations per frame
    float remap = 0.00007;
    // if it takes e.g. 10 seconds for 1 year, then 600 frames = TIMELIM/DT iters = 24*365 = 8760
    // 8760/600 = 14 iterations per frame
    while (!WindowShouldClose()) {
        while (t < TIMELIM) {
            for (int subiter = 0; subiter < subLim; subiter++) {
                for (int i = 0; i < BODYCNT; i++) {
                    // f = g m1 m2 / r^2
                    // a = f/m
                    Vector3 acc = Vector3Zero();
                    for (int j = 0; j < BODYCNT; j++) {
                        if (i == j)
                            continue;
                        double r2 = Vector3DistanceSqr(planets[i].pos, planets[j].pos);
                        // int f = G * planets[i].mass * planets[j].mass / (r2);
                        float a = G * planets[j].mass / (r2);
                        if (a == 0.0 && iter % 365 == 0) {
                            printf("i:%d j:%d WARNING: a is 0 \n", i, j);
                        }
                        acc = Vector3Add(acc,
                                         Vector3Scale(
                                             Vector3Normalize(
                                                 Vector3Subtract(
                                                     planets[j].pos, planets[i].pos)),
                                             a));
                    }
                    // f=ma, a=f/m
                    planets[i].acc = acc;
                    planets[i].vel = Vector3Add(planets[i].vel, Vector3Scale(planets[i].acc, DT));
                    planets[i].pos = Vector3Add(planets[i].pos, Vector3Scale(planets[i].vel, DT));

                    if (Vector3LengthSqr(planets[i].vel) == 0 && i != 0) {
                        puts("warning: velocity is 0");
                    }

                    if ((iter % (50000 * 60)) == 0) {
                        if (i == 4) {
                            printf("iter: %ld jupiter x %f y %f\n", iter, planets[i].pos.x, planets[i].pos.y);
                        }

                        printf("planet %d %s\n", i, planets[i].name);
                        printf("distance from sun: %f AU\n", Vector3Length(planets[i].pos) / AU);
                        if (i == 3) {
                            printf("moon distance from earth: %f AU\n", Vector3Length(Vector3Subtract(planets[i].pos, planets[2].pos)) / AU);
                        }
                        // printf("position in AU: (%f, %f, %f)\n", planets[i].pos.x / AU, planets[i].pos.y / AU, planets[i].pos.z / AU);
                        // printf("velocity: (%f, %f, %f)\n", planets[i].vel.x, planets[i].vel.y, planets[i].vel.z);
                        // printf("acceleration: (%f, %f, %f)\n", planets[i].acc.x, planets[i].acc.y, planets[i].acc.z);
                    }
                }
                iter++;
            }

            BeginDrawing();
            ClearBackground(RAYWHITE);
            // find max x and y to scale
            float maxX = 0.;
            float maxY = 0.;
            for (int i = 0; i < BODYCNT; i++) {
                if (abs(planets[i].pos.x) > maxX)
                    maxX = planets[i].pos.x;
                if (abs(planets[i].pos.y) > maxY)
                    maxY = planets[i].pos.y;
            }
            // scale to fit screen given fixed aspect ratio of say 800/450 = 1.78
            const float padding = 100.0;
            maxXY = max(maxY * (800.f / 450.f), maxX);
            // remap = max(.0000000005, min(remap, .5 * ((800.0 - padding) / maxXY)));
            remap = max(.00000005, min(remap, .5 * ((800.0 - padding) / (maxXY * .05))));

            if (GetMouseWheelMove() != 0) {
                printf("mouse wheel move %f, remap was %f\n", GetMouseWheelMove(), remap);
                if (GetMouseWheelMove() < 0.) {
                    remap *= (1 / (-GetMouseWheelMove() + .5));
                    printf("remap is now %f\n", remap);
                } else {
                    remap *= min(3., GetMouseWheelMove() + .5);
                }
            }
            // printf("scale: %f,maxXY %f \n", scale, maxXY);

            for (int i = 0; i < BODYCNT; i++) {
                // maybe draw trails?
                Planet activePlanet = planets[2];
                Planet planet = planets[i];
                Vector3 transPlanetPos = Vector3Subtract(
                    Vector3Add(
                        Vector3Scale(planet.pos,
                                     remap),
                        centre),
                    Vector3Scale(activePlanet.pos,
                                 remap)); // can change frame of reference
                Vector2 flatPos = {transPlanetPos.x, transPlanetPos.y};
                // printf("planet %s x %f y %f radius %f \n", planet.name, flatPos.x, flatPos.y, planet.radius * (scale * 1000));
                DrawCircleV(flatPos, planet.radius * (remap * 1000), planet.color);
                // can use drawcirclegradient for 2 cols!
                // not to scale radius because they would be tiny
            }

            DrawFPS(10, 10);
            EndDrawing();

            t += DT;
        }
        puts("done");
        CloseWindow();
        return 0;
    }
    CloseWindow();
    return 0;
}