#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "rlFPCamera/rlFPCamera.h"

#define RAYGUI_IMPLEMENTATION

#include "raygui.h"

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define clamp(a, b, c) (max(b, min(a, c)))

void UpdateDrawFrame(void);

#if defined(PLATFORM_WEB)
#include <emscripten/emscripten.h>
#define TRAIL_LEN 1e5
#else
#define TRAIL_LEN 1e7
#endif
#define XYSIZE .5 // remember this has quadratic memory usage!
// I build with cc nbody_chart.c -O3 -Wall -Wextra -pedantic -ffast-math
// -funsafe-math-optimizations -l:libraylib.a -lm -pthread -o nbody_chart
// from
// https://stackoverflow.com/questions/2509679/how-to-generate-a-random-integer-number-from-within-a-range
unsigned int rand_interval(unsigned int min, unsigned int max) {
  int r;
  const unsigned int range = 1 + max - min;
  const unsigned int buckets = RAND_MAX / range;
  const unsigned int limit = buckets * range;

  /* Create equal size buckets all in a row, then fire randomly towards
   * the buckets until you land in one of them. All buckets are equally
   * likely. If you land off the end of the line of buckets, try again. */
  do {
    r = rand();
  } while (r >= limit);

  return min + (r / buckets);
}

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
// challenge 2: plot elliptical orbits of planets (Done)
// assume sun is stationary origin, orbits are ellipses in same plane
// challenge 3: create a 2d animation of the planets orbiting the sun (Done)
// plot inner 5 planets seperately from outer planets
// inner planets have 1 earth orbit per s, outer planets have 1 jupiter orbit
// per s
// challenge 4: create a 3d animation of the planets orbiting the sun
// Pluto is basically the only one which will be noticeably inclined, but add
// data for all! (done)
// challenge 5: calculate orbit angle vs time for an eccentric orbit
// e.g. Pluto and compare to a circular one with same period will have to
// do numerical integration with cumsum, simpson's rule
// challenge 6: solar system spirograph
// Choose a pair of planets and determine their orbits vs time.
// At time intervals of Dt, draw a line between the planets and plot this
// line. Keep going for N orbits of the outermost planet.
// challenge 7: like Ptolemy, use orbital models to plot the orbits of bodies relative to another
// orbiting body, e.g. earth
// extensions! make n body sims of binary stars! do
// more funky exoplanet stuff data can be found at visualise lagrange points
// relativistic corrections
// {\displaystyle \sigma ={\frac {24\pi
// ^{3}L^{2}}{T^{2}c^{2}\left(1-e^{2}\right)}}\ ,} sigma = (24 pi^3
// a^2)/(T^2c^2(1-e^2)) where sigma is the perihelion shift in radians per
// revolution (only really relevant for Mercury and Venus)

// https://www.youtube.com/watch?v=Qw6uZQYwG7A
// write up a report/paper in latex
// can do wasm stuff to make a website!
// visualise transfer orbits (hohmann, const thrust (nuclear, solar sails, ion
// engines)), gravity assists, etc. variables e (normally \epsilon) =
// eccentricity, a = semi-major axis, b = semi-minor axis,T = orbital period,
// beta = orbital inclination in degrees

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
// because (M_sun + m_planet) ≈ M_sun, we can ignore mass of planet, finding the
// following equation: (orbital period/Earth Year) ≈ (a/AU)^3/2 (This can be
// graphed for both our solar system and exoplanet.eu data) line of best fit can
// be found and grad diff from linear should be v small (remember, best fit
// should go through 0,0)
// https://en.wikipedia.org/wiki/List_of_multiplanetary_systems

// there are various ways to detect exoplanets (see here:
// https://en.wikipedia.org/wiki/Methods_of_detecting_exoplanets)

// semi major axis can be found by max dist from barycentre
// period can be found by time to start pos within epsilon e.g. (in AU) 1e-6
// for data:

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
  RenderTexture2D label;
  Model object;
  Color color; // maybe later add texture, atmosphere, etc.
} Body;

typedef Body Planet; // "OOP"
typedef Body Moon;
typedef Body Star;
const float scale = 1e-3; // check this
const float radScale = 1;
const float AU = 1.496e11 * scale; // delta of 1e4
float dt = 10000.;                 // DT of 10000 = 3hrs per frame, ~week per second
float currTime = 0;
unsigned int frameCnt = 0;
const size_t subLim = 1; // steps per frame
float remap = 1e-8;
RenderTexture2D target;
Model plane;
bool XorZ = false;
bool spiro = false;
int spiroStep = 2;
bool help = true;
bool axes = true;
bool paused = false;
bool trails = false;
bool cursorState = false;
Vector3 StartPosition;
float fovy;
rlFPCamera camera;
Camera3D *ViewCamera;
// size_t iter = 0;
Star sun = {
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
    .e = 0.20563069,
    .inclination = 7. * PI / 180.,
    .radius = 6.371e6,
    .color = LIGHTGRAY,
    .parent = &sun,
};
Planet venus = {
    .name = "Venus",
    .a = 0.72333199 * AU,
    .e = 0,
    .inclination = 3.4 * PI / 180.,
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
Moon moon = {
    .name = "Moon",
    .a = 3.844e8 * scale * 50, // not to scale; hand tuned because earth radius was inflated
    .e = 0,
    .inclination = 5.1 * PI / 180.,
    .radius = 1.7374e6,
    .color = GRAY,
    .parent = &earth,
};
Planet mars = {
    .name = "Mars",
    .a = 1.52366231 * AU,
    .e = 0,
    .inclination = 1.85 * PI / 180.,
    .radius = 3.3895e6,
    .color = RED,
    .parent = &sun,
};
Moon phobos = {
    .name = "Phobos",
    .a = 9.378e7 * scale * 50, // not to scale; hand tuned because planetary radius was inflated
    .e = 0,
    .inclination = 0,
    .radius = 1.1267e4 * 10,
    .color = GRAY,
    .parent = &mars,
};
Moon deimos = {
    .name = "Deimos",
    .a = 2.3e8 * scale * 50, // not to scale; hand tuned because planetary radius was inflated
    .e = 0,
    .inclination = 0,
    .radius = 6.2e3 * 10,
    .color = GRAY,
    .parent = &mars,
};
Planet jupiter = {
    .name = "Jupiter",
    .a = 5.204267 * AU,
    .e = 0,
    .inclination = 1.3 * PI / 180.,
    .radius = 6.9911e7,
    .color = ORANGE,
    .parent = &sun,
};
Planet saturn = {
    .name = "Saturn",
    .a = 9.582017 * AU,
    .e = 0,
    .inclination = 2.5 * PI / 180.,
    .radius = 5.8232e7,
    .color = BROWN,
    .parent = &sun,

};
Planet uranus = {
    .name = "Uranus",
    .a = 19.22941195 * AU,
    .e = 0,
    .inclination = 0.8 * PI / 180.,
    .radius = 2.5362e7,
    .color = SKYBLUE,
    .parent = &sun,
};
Planet neptune = {
    .name = "Neptune",
    .a = 30.10366151 * AU,
    .e = 0,
    .inclination = 1.8 * PI / 180.,
    .radius = 2.4622e7,
    .color = DARKBLUE,
    .parent = &sun,
};
Planet pluto = {
    .name = "Pluto",
    .a = 39.48211675 * AU,
    .e = .244,
    .inclination = 17 * PI / 180.,
    .radius = 1.188e6,
    .color = DARKGRAY,
    .parent = &sun,
};
Body *bodies[] = {&sun,    &mercury, &venus,  &earth,  &moon,    &mars, &phobos,
                  &deimos, &jupiter, &saturn, &uranus, &neptune, &pluto};
size_t bodycnt = sizeof(bodies) / sizeof(bodies[0]);

int centrePlanetDropdownActive = 0;
bool centrePlanetDropdownEditMode = false;

Body *firstSpiroPlanet;
int firstSpiroDropdownActive = 0;
bool firstSpiroDropdownEditMode = false;

Body *secondSpiroPlanet;
int secondSpiroDropdownActive = 0;
bool secondSpiroDropdownEditMode = false;


Planet *CentrePlanet = &sun;
Vector3 starPos[100] = {0};
#ifdef PLATFORM_WEB
  int width = 1200;
  int height = 500;
#else
  int width = 1900;
  int height = 900;
#endif
int main(void) {
  SetConfigFlags(FLAG_MSAA_4X_HINT);
  InitWindow(width, height, "raylib planets - basic demo");
  // skybox
  // https://github.com/petrocket/spacescape
  // Mesh cube = GenMeshCube(1.0f, 1.0f, 1.0f);
  // Model skybox = LoadModelFromMesh(cube);
  // Shader sky_shader = LoadShader(0, "skybox_shader.glsl");
  // skybox.materials[0].shader = sky_shader;

  // SetShaderValue(skybox.materials[0].shader,
  // GetShaderLocation(skybox.materials[0].shader, "environmentMap"), (int[1]){
  // MATERIAL_MAP_CUBEMAP }, SHADER_UNIFORM_INT);
  // SetShaderValue(skybox.materials[0].shader,
  // GetShaderLocation(skybox.materials[0].shader, "doGamma"), (int[1]) { useHDR
  // ? 1 : 0 }, SHADER_UNIFORM_INT); SetShaderValue(skybox.materials[0].shader,
  // GetShaderLocation(skybox.materials[0].shader, "vflipped"), (int[1]){ useHDR
  // ? 1 : 0 }, SHADER_UNIFORM_INT);

  plane = LoadModelFromMesh(GenMeshPlane(XYSIZE * 10.0f, XYSIZE * 10.0f, 1, 1));
  unsigned int minStarDist = 40 * AU * remap;
  for (int i = 0; i < 100; i++) {
    signed int randSign = ((rand() % 2 == 0) ? -1 : 1);
    starPos[i] = (Vector3){
        (float)((rand() % 2 == 0) ? -1 : 1) * rand_interval(minStarDist, 1.5 * minStarDist),
        (float)((rand() % 2 == 0) ? -1 : 1) * rand_interval(minStarDist, 1.5 * minStarDist),
        (float)((rand() % 2 == 0) ? -1 : 1) * rand_interval(minStarDist, 1.5 * minStarDist)};
  }
  Vector3 StartPosition = (Vector3){1.0f, 0.0f, 0.0f}; // Camera position
  float fovy = 45.0f;                                  // Camera field-of-view Y
  rlFPCameraInit(&camera, fovy, StartPosition);
  ViewCamera = &(camera.ViewCamera);
  ViewCamera->up = (Vector3){0.0f, 1.0f, 0.0f};
  ViewCamera->target = Vector3Zero();
  camera.MoveSpeed.z = 10;
  camera.MoveSpeed.x = 5;
  camera.FarPlane = 5000;
  camera.AllowFlight = true;
  for (size_t i = 0; i < bodycnt; i++) {
    Planet *body = bodies[i];
    // allocate memory for trail
    body->trail.data = malloc(sizeof(Vector3) * TRAIL_LEN);
    if (body->trail.data == NULL) {
      fputs("Out of Memory!\n", stderr);
      return 1; // OOM, OS can do cleanup i guess
    }
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
    else {
      body->b = sqrtf(body->a * body->a * (1 - body->e * body->e));
      printf("%s :a:%f, b:%f\n", body->name, body->a, body->b);
    }
    body->label = LoadRenderTexture(400 * XYSIZE, 800 * XYSIZE);
    BeginTextureMode(body->label);
    ClearBackground(BLANK);
    DrawText(body->name, 0, 0, 100 * XYSIZE, BLACK);
    EndTextureMode();
  }

#if defined(PLATFORM_WEB)
  emscripten_set_main_loop(UpdateDrawFrame, 0, 1);
#else
  SetTargetFPS(60);
  while (!WindowShouldClose()) // Detect window close button or ESC key
  {
    UpdateDrawFrame();
  }
#endif
  // cleanup:
  for (size_t i = 0; i < bodycnt; i++) {
    Planet *body = bodies[i];
    free(body->trail.data);
    UnloadRenderTexture(body->label);
  }
}
void UpdateDrawFrame(void) {
  BeginDrawing();

  // TODO: add starfield option
  // https://www.shadertoy.com/view/Md2SR3
  // https://github.com/raysan5/raylib/blob/bc9c06325481c0b4b5a2db3b2a8281465569ba3e/examples/models/models_skybox.c

  ClearBackground(RAYWHITE);
  // If you change ClearBackground col, make sure to change the text color
  if (!cursorState){
  rlFPCameraUpdate(&camera);
}
  if (GetMouseWheelMove() != 0) {
    ViewCamera->fovy += GetMouseWheelMove();
    ViewCamera->fovy = clamp(355, 1, ViewCamera->fovy);
  }
  if (IsKeyPressed(KEY_R)) {
    camera.CameraPosition = StartPosition;
    ViewCamera->target = Vector3Zero();
    // ViewCamera->fovy = fovy;
    dt = 10000.;
    CentrePlanet = &sun;
  }
  if (IsKeyPressed(KEY_SPACE)) {
    paused = !paused;
  }
  if (IsKeyPressed(KEY_X)) {
    axes = !axes;
  }
  if (IsKeyPressed(KEY_T)) {
    trails = !trails;
  }
  if (IsKeyPressed(KEY_COMMA)) {
    dt *= (1.0 / 1.1);
  }
  if (IsKeyPressed(KEY_PERIOD)) {
    dt *= 1.1;
  }
  if (IsKeyPressed(KEY_H)) {
    help = !help;
  }
  if (IsKeyPressed(KEY_X)) {
    XorZ = !XorZ;
  }
  if (IsKeyPressed(KEY_P)) {
    CentrePlanet = &earth;
  }
  if (IsKeyPressed(KEY_L)) {
    spiro = !spiro;
  }
  if (IsKeyPressed(KEY_O)) {
    // topdown
    camera.CameraPosition = (Vector3){0, 20, 0};
    ViewCamera->target = (Vector3){0, 19, 0};
  }
  if (IsKeyPressed(KEY_K)){
    if (!cursorState){
    ShowCursor();
    EnableCursor();
    }
    else
    {
    HideCursor();
    DisableCursor();
    }
    cursorState = !cursorState;
  }
  const char* bodylist = "Sun;Mercury;Venus;Earth;Moon;Mars;Jupiter;Saturn;Uranus;Neptune;Pluto";
  GuiUnlock();
  GuiSetStyle(DROPDOWNBOX,TEXT_ALIGNMENT,TEXT_ALIGN_CENTER);
  DrawText("Centre Planet:",width-200,20,20,BLACK);
  if (GuiDropdownBox((Rectangle){width-200,40,125,30},bodylist,&centrePlanetDropdownActive,centrePlanetDropdownEditMode))
  {
    centrePlanetDropdownEditMode = !centrePlanetDropdownEditMode;
    CentrePlanet = bodies[centrePlanetDropdownActive]; 
  }
  DrawText("Spirograph Planets:",width-500,20,20,BLACK);
  if (GuiDropdownBox((Rectangle){width-500,40,125,30},bodylist,&firstSpiroDropdownActive,firstSpiroDropdownEditMode)){
    firstSpiroDropdownEditMode = !firstSpiroDropdownEditMode;
    firstSpiroPlanet = bodies[firstSpiroDropdownActive];
  }
  if (GuiDropdownBox((Rectangle){width-500,80,125,30},bodylist,&secondSpiroDropdownActive,secondSpiroDropdownEditMode)){
    secondSpiroDropdownEditMode = !secondSpiroDropdownEditMode;
    secondSpiroPlanet = bodies[secondSpiroDropdownActive];
  }
  GuiLock();

if (!paused){
    frameCnt++;
  }
  for (size_t subiter = 0; subiter < subLim * (size_t)(dt / 10000.); subiter++) {
    if (subiter == 0) {
      // prevent drawing once per Body by drawing outside of body loop
      // draw iteration number in days (ITER*DT/3600/24)
      DrawText(TextFormat("Days: %.1f", currTime / 3600. / 24.), 10, 30, 20, BLACK);
      DrawText(TextFormat("Years: %.1f", currTime / 3600. / 24. / 365.), 10, 50, 20, BLACK);
      // Draw DT in hours
      DrawText(TextFormat("Hours per second: %.1f", dt / 3600. * 60.), 10, 70, 20, BLACK);
      DrawText(TextFormat("FOV: %.1f", ViewCamera->fovy), 10, 90, 20, BLACK);

      Vector3 pos = ViewCamera->position;
      Vector3 target = ViewCamera->target;
      DrawText(TextFormat("pos: x:%.1f y:%.1f z:5.1f", pos.x, pos.y, pos.z), 10, 110, 20, BLACK);
      // DrawText(TextFormat("target: x:%.1f y:%.1f
      // z:5.1f",target.x,target.y,target.z),10,110,20,BLACK);
      if (help) {
        DrawText(
            "Help info: press L to toggle axes, H to toggle help, \n[Space] to toggle pause, T "
            "to toggle trails\n[comma] and [period] to slow down and speed up the sim respectively",
            10, 130, 20, BLACK);
      }
    }
    rlFPCameraBeginMode3D(&camera);
    if (!paused) {
      currTime += dt;
    }
    if (axes) {
      DrawGrid(10, 1.0f);
      DrawLine3D(Vector3Zero(), (Vector3){0, 0, 1000.f}, BLUE);
      DrawLine3D(Vector3Zero(), (Vector3){0, 1000.f, 0}, BLUE);
      DrawLine3D(Vector3Zero(), (Vector3){1000.f, 0, 0}, BLUE);
    }
    for (int i = 0; i < 100; i++) {
      // DrawSphere(starPos[i], 1, WHITE);
    }
    for (size_t i = 0; i < bodycnt; i++) {
      Planet *body = bodies[i];

      Vector3 pos = {
          body->a * cos(body->theta) * cos(body->inclination),
          body->a * sin(body->theta) * sin(body->inclination),
          body->b * sin(body->theta),
      };

      if (body->parent != NULL) {
        pos = Vector3Add(pos, body->parent->pos);
      }
      body->pos = pos;
      Vector3 OffsetPos = Vector3Subtract(pos, CentrePlanet->pos);
      Vector3 mappedPos = Vector3Scale(OffsetPos, remap);

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
        // rotate to match orbital plane
        plane.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = body->label.texture;
        rlDisableBackfaceCulling();
        // can calculate rotation with dot product
        Vector3 diffvec = Vector3Subtract(ViewCamera->position, mappedPos);
        Vector3 xydiff = Vector3Multiply((Vector3){1, 1, 0}, diffvec);
        float xAngle = Vector3Angle(xydiff, (Vector3){1.0f, 0.0f, 0.0f});
        float zAngle = Vector3Angle(xydiff, (Vector3){0.0f, 1.0f, 0.0f});
        Vector3 scale = Vector3Scale(Vector3One(), 5e-2);
        // printf("%s:\t%f %f %f\n", body->name, xAngle, yAngle, xAngle - yAngle);
        // printf("xy offsets:\t%f %f\n", camera.position.x, camera.position.y);
        float angle = 0;
        // I cba to implement labels, too fiddly
        /* rlRotatef(90,1,0,0); // rotate plane to be vertical not flat on floor
         DrawModelEx(plane, Vector3Add(mappedPos, (Vector3){0, body->radius * remap * 1.5f, 0}),
                     (Vector3){1.0f, 0.0f, 0.0f}, angle * RAD2DEG, scale, WHITE);
        */
        // tint is white because tint is multiplied with diffuse color
        rlEnableBackfaceCulling();

        DrawSphere(mappedPos, body->radius * remap, body->color);
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
            Vector3 pos1 = body->trail.data[j];
            Vector3 pos2 = body->trail.data[(j + 1) % len];
            // Vector3 OffsetPos1 = Vector3Subtract(pos1,CentrePlanet->pos);
            // Vector3 OffsetPos2 = Vector3Subtract(pos2,CentrePlanet->pos);
            Vector3 OffsetPos1 = Vector3Subtract(pos1, CentrePlanet->trail.data[j]);
            Vector3 OffsetPos2 = Vector3Subtract(pos2, CentrePlanet->trail.data[(j + 1) % len]);
            Vector3 mappedPos1 = Vector3Scale(OffsetPos1, remap);
            Vector3 mappedPos2 = Vector3Scale(OffsetPos2, remap);
            DrawLine3D(mappedPos1, mappedPos2, body->color);
            if (spiro && body == firstSpiroPlanet && (j % spiroStep == 0)) {
              Vector3 otherPos = secondSpiroPlanet->trail.data[j];
              Vector3 offsetOtherPos = Vector3Subtract(otherPos, CentrePlanet->trail.data[j]);
              Vector3 mappedOtherPos = Vector3Scale(offsetOtherPos, remap);
              DrawLine3D(mappedPos1, mappedOtherPos, GREEN);
            }
          }
        }
      }
      // t ^ 2 = k *a ^ 3
      // k=t^2/a^3 = (365*3600*24)^2/(AU)^3 = 2.9704e-19
      double k = 2.9704e-19 * 1 / (scale * scale * scale); // scale down AU
      double T = sqrt(k * body->a * body->a * body->a);
      if (frameCnt == 10) {
        printf("body:%s has T (in years) of:%1f\n", body->name, T / 3600. / 24. / 365.);
      }
      if (T != 0 && !paused) { // don't divide by 0 (Sun)
        // body->theta += dt / T / (subLim * (size_t)(dt / 10000.)) * 2 * PI;
        body->theta += (dt / T) / subLim * 2 * PI;
      }
    }
    // draw spirograph between 2 planets

    rlFPCameraEndMode3D();
  }
  DrawFPS(10, 10);
  EndDrawing();
}
