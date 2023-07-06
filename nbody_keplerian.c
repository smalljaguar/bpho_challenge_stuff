#include "raylib.h"
#include "raymath.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define clamp(a, b, c) (max(b, min(a, c)))
#define TRAIL_LEN 1e8
// I build with cc nbody_chart.c -O3 -Wall -Wextra -pedantic -ffast-math -funsafe-math-optimizations -l:libraylib.a -lm -pthread -o nbody_chart

// more "deterministic" than sim
// use parametric equation for ellipse given semi major and semi minor axis a,b
// https://math.stackexchange.com/questions/3127984/whats-the-parametric-equation-of-a-partial-ellipse-in-3d-space-with-given-major
// z = c sin(theta)
// y = b sin(theta)
// x = a cos(theta)

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

typedef struct CircularArray {
    size_t len;
    size_t pos;
    size_t start;
    Vector3 *data;
} CircularArray;

typedef struct Body {
    // probably should use SI units I guess
    float a;             // semi major axis
    float b;             // semi minor axis
    float e;             // eccentricity
    float inclination;   // inclination in degrees/radians?
    CircularArray trail; // array of vec3s for trail
    Vector3 pos;
    double radius;
    double theta;
    const char *name;
    struct Body *parent;
    Color color; // maybe later add texture, atmosphere, etc.
} Body;

typedef Body Planet; // "OOP"

#define BODYCNT 11
int main() {
    const float scale = 1e-3;
    const float radScale = 1;
    const float AU = 1.496e11 * scale;              // delta of 1e4
    float dt = 10000.;                              // DT of 10000 = 3hrs per frame, ~week per second
    const float TIMELIM = 3600. * 24. * 365 * 100.; // 100 years
    const size_t subLim = 1;                        // steps per frame
    float remap = 5e-6;
    float mult = 0;
    int labels = 1;
    bool paused = false;
    bool trails = false;
    // size_t iter = 0;
    Body sun = {
        .name = "Sun",
        .a = 0,
        .e = 0,
        .inclination = 0,
        .radius = 6.957e8 * 5e-2, // sun is too big!!
        .color = RED,
        .parent = NULL, // root node
    };
    Planet mercury = {
        .name = "Mercury",
        .a = 0.38709893 * AU,
        .e = 0,
        .inclination = 0,
        .radius = 6.371e6,
        .color = LIGHTGRAY,
        .parent = &sun,
    };
    Planet venus = {
        .name = "Venus",
        .a = 0.72333199 * AU,
        .e = 0,
        .inclination = 0,
        .radius = 6.051e6,
        .color = YELLOW,
        .parent = &sun,
    };

    Planet earth = {
        .name = "Earth",
        .a = AU,
        .e = 0,
        .inclination = 0,
        .radius = 6.371e6,
        .color = BLUE,
        .parent = &sun,
    };
    Body moon = {
        .name = "Moon",
        .a = 3.844e8 * scale * 50, // not to scale; hand tuned because earth radius was inflated
        .e = 0,
        .radius = 1.7374e6,
        .color = GRAY,
        .parent = &earth,
    };
    Planet mars = {
        .name = "Mars",
        .a = 1.52366231 * AU,
        .e = 0,
        .inclination = 0,
        .radius = 3.3895e6,
        .color = RED,
        .parent = &sun,
    };
    Body phobos = {
        .name = "Phobos",
        .a = 9.378e7 * scale * 50, // not to scale; hand tuned because planetary radius was inflated
        .e = 0,
        .radius = 1.1267e4 * 10,
        .color = GRAY,
        .parent = &mars,
    };
    Body deimos = {
        .name = "Deimos",
        .a = 2.3e8 * scale * 50, // not to scale; hand tuned because planetary radius was inflated
        .e = 0,
        .radius = 6.2e3 * 10,
        .color = GRAY,
        .parent = &mars,
    };
    Planet jupiter = {
        .name = "Jupiter",
        .a = 5.204267 * AU,
        .e = 0,
        .inclination = 0,
        .radius = 6.9911e7,
        .color = ORANGE,
        .parent = &sun,
    };
    Planet saturn = {
        .name = "Saturn",
        .a = 9.582017 * AU,
        .e = 0,
        .inclination = 0,
        .radius = 5.8232e7,
        .color = BROWN,
        .parent = &sun,

    };
    Planet pluto = {
        .name = "Pluto",
        .a = 39.48211675 * AU,
        .e = 0,
        .inclination = 0,
        .radius = 1.188e6,
        .color = DARKGRAY,
        .parent = &sun,
    };
    InitWindow(800, 450, "raylib planets - basic demo?");
    Vector2 offset = {(float)GetScreenWidth() / 2, (float)GetScreenHeight() / 2};
    SetTargetFPS(60);

    Body *bodies[BODYCNT] = {&sun, &mercury, &venus, &earth, &moon, &mars,
                             &phobos, &deimos, &jupiter, &saturn, &pluto};

    for (int i = 0; i < BODYCNT; i++) {
        Planet *body = bodies[i];
        // allocate memory for trail
        body->trail.data = malloc(sizeof(Vector3) * TRAIL_LEN);
        body->trail.len = TRAIL_LEN;
        body->trail.start = 0;
        body->trail.pos = 0;
        // scale down Radii at runtime, calculate b to prevent error
        body->radius = body->radius * radScale;
        // e = sqrt(1-(b^2/a^2))
        // e^2*a^2 = a^2 - b^2
        // b^2 = a^2(1-e^2)
        if (body->a == 0 || body->e == 1)
            body->b = 0;
        else
            body->b = sqrt(body->a * body->a * (1 - body->e * body->e));
    }
    for (float currTime = 0; currTime < TIMELIM; currTime += dt) {
        BeginDrawing();
        ClearBackground(RAYWHITE);
        if (GetMouseWheelMove() != 0) {
            printf("mouse wheel move %f, remap was %f\n", GetMouseWheelMove(), remap);
            if (GetMouseWheelMove() < 0.) {
                mult = max(.2, (1 / (-GetMouseWheelMove() + .5)));
            } else {
                mult = min(5., GetMouseWheelMove() + .5);
            }
            remap *= mult;
            printf("remap is now %f\n", remap);
        }
        if (IsKeyDown(KEY_LEFT)) {
            offset.x += 10;
        }
        if (IsKeyDown(KEY_RIGHT)) {
            offset.x -= 10;
        }
        if (IsKeyDown(KEY_UP)) {
            offset.y += 10;
        }
        if (IsKeyDown(KEY_DOWN)) {
            offset.y -= 10;
        }
        if (IsKeyPressed(KEY_R)) {
            offset = (Vector2){(float)GetScreenWidth() / 2, (float)GetScreenHeight() / 2};
            remap = 5e-6;
            dt = 10000.;
        }
        if (IsKeyPressed(KEY_SPACE)) {
            paused = !paused;
        }
        if (IsKeyPressed(KEY_L)) {
            labels = !labels;
        }
        if (IsKeyPressed(KEY_T)) {
            trails = !trails;
        }
        if (IsKeyPressed(KEY_COMMA)) {
            dt *= .9;
        }
        if (IsKeyPressed(KEY_PERIOD)) {
            dt *= 1.1;
        }
        if (IsKeyPressed(KEY_ESCAPE)) {
            goto end;
        }
        for (size_t subiter = 0; subiter < subLim * (size_t)(dt / 10000.); subiter++) {
            for (int i = 0; i < BODYCNT; i++) {
                Planet *body = bodies[i];

                // Vector3 pos = {
                //     body.a * cos(body.theta) * cos(body.inclination),
                //     body.b * sin(body.theta),
                //     body.a * sin(body.theta) * cos(body.inclination),
                // };
                Vector3 pos = {
                    body->a * cos(body->theta),
                    body->b * sin(body->theta),
                    0.,
                };
                if (body->parent != NULL) {
                    pos = Vector3Add(pos, body->parent->pos);
                }
                body->pos = pos;
                Vector3 offset3 = {offset.x, offset.y, 0};
                Vector3 mappedPos = Vector3Add(
                    Vector3Scale(body->pos,
                                 remap),
                    offset3);
                Vector2 flatPos = {mappedPos.x, mappedPos.y};

                // add to trail (should it be every frame or subiter?)
                if (!paused) {
                    body->trail.data[body->trail.pos] = body->pos;

                    // circular!
                    if (body->trail.pos >= body->trail.len - 1) {
                        body->trail.start++;
                        body->trail.pos = 0;
                    } else {
                        body->trail.pos++;
                    }
                    if (body->trail.start != 0) {
                        // overwriting old trails
                        body->trail.start++;
                    }
                    if (body->trail.start >= body->trail.len - 1) {
                        body->trail.start = 0;
                    }
                }
                if (subiter == 0) {
                    // only draw once per frame
                    if (labels) {
                        DrawText(body->name, flatPos.x,
                                 flatPos.y + body->radius * remap, 10, BLACK);
                    }
                    DrawCircleV(flatPos, body->radius * (remap), body->color);
                    // draw iteration number in days (ITER*DT/3600/24)
                    DrawText(TextFormat("Days: %.1f", currTime / 3600. / 24.),
                             10, 30, 20, BLACK);
                    // Draw DT in hours
                    DrawText(TextFormat("Hours per second: %.1f", dt / 3600. * 60.), 10, 50, 20, BLACK);
                    // draw trail
                    // TODO: fix drawing after circular overwrite 
                    if (trails) {
                        size_t start = body->trail.start;
                        size_t end = body->trail.pos - 1; // to account for pos2 being j+1
                        size_t len = body->trail.len;
                        size_t usedLen = end - start;
                        // probably festooned with off-by-one errors
                        if (end > start) {
                            usedLen = end - start;
                        } else {
                            usedLen = len - 1;
                        }
                        for (size_t j = 0; j < usedLen; j++) {
                            j = (j + start) % len;
                            Vector2 pos1 = {body->trail.data[j].x, body->trail.data[j].y};
                            Vector2 pos2 = {body->trail.data[(j + 1) % len].x, body->trail.data[(j + 1) % len].y};
                            pos1 = Vector2Add(Vector2Scale(pos1, remap), offset);
                            pos2 = Vector2Add(Vector2Scale(pos2, remap), offset);
                            DrawLineV(pos1, pos2, body->color);
                        }
                    }
                }
                // t ^ 2 = k *a ^ 3
                // k=t^2/a^3 = (365*3600*24)^2/(AU)^3 = 2.9704e-19
                double k = 2.9704e-19 * 1 / (scale * scale * scale); // scale down AU
                double T = sqrt(k * body->a * body->a * body->a);
                if (T != 0 && !paused) { // don't divide by 0 (Sun)
                    body->theta += dt / T / (subLim * (size_t)(dt / 10000.)) * 2 * PI;
                }
            }
        }
        DrawFPS(10, 10);
        EndDrawing();
        if (WindowShouldClose())
            goto end; // this is done to exit current block scope and so avoid GLFW errors
        // Don't need to free because basically everything is on the stack
        // kinda jank but works
    }
end:
    CloseWindow();
}