#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdbool.h>

#include <fftw3.h>
#include <math.h>
#include <SDL/SDL.h>

#define PALETTE_LEN_BITS 12
#define PALETTE_LEN (1 << PALETTE_LEN_BITS)

typedef struct {
  Uint32 *colors;
  unsigned int len;
} palette_t;

typedef struct {
  float pos;
  float r;
  float g;
  float b;
} palette_point_t;

static void *malloc_check(size_t len) {
  void *p;
  p = malloc(len);
  if (! p) {
    printf("No mem.\n");
    exit(-1);
  }
  return p;
}

void set_color(palette_t *palette, int i, float r, float g, float b,
               SDL_PixelFormat *format) {
  if (i >= palette->len)
    return;
  palette->colors[i] = SDL_MapRGB(format, r * 255, g * 255, b * 255);
}

/* Generates a color palette, setting palette->colors and palette->len.
 * Allocates new memory for palette->colors (is not freed or reallocd).
 * 'n_colors' defines how many colors are generated in the palette.
 * 'points' is a definition colors at specific intervals, 'n_points' gives the
 * number of palette_point_t array elements in 'points'.
 * 'format' is used to generate video mode specific color data. */
void make_palette(palette_t *palette, int n_colors,
                  palette_point_t *points, int n_points,
                  SDL_PixelFormat *format) {
  int i;

  palette->colors = malloc_check(n_colors * sizeof(Uint32));
  palette->len = n_colors;


  if (n_points < 1) {
    for (i = 0; i < palette->len; i++) {
      float val = (float)i / palette->len;
      set_color(palette, i, val, val, val, format);
    }
    return;
  }

  palette_point_t *last_p = points;
  palette_point_t *first_p = points;

  for (i = 1; i < n_points; i ++) {
    if (points[i].pos > last_p->pos)
      last_p = &points[i];
    if (points[i].pos < first_p->pos)
      first_p = &points[i];
  }
  if (last_p->pos > 1.0) {
    float norm_factor = last_p->pos;
    for (i = 0; i < n_points; i ++)
      points[i].pos /= norm_factor;
  }
  
  // duplicate the last point to "the left", wrap back below zero.
  palette_point_t p = *last_p;
  p.pos -= 1.0;
  // ...unless another point is defined there.
  if (p.pos >= first_p->pos)
    p = *first_p;

  // also duplicate the first point to "the right".
  palette_point_t post_last = *first_p;
  post_last.pos += 1.0;

  int color_pos = 0;

  while(color_pos < n_colors) {

    // look for the next point, the one with the next largest pos after p.pos
    palette_point_t *next_p = NULL;
    for (i = 0; i < n_points; i ++) {
      float i_pos = points[i].pos;
      if ((i_pos > p.pos)
          &&
          (
           (! next_p)
           || (i_pos < next_p->pos)
          )
         )
        next_p = &points[i];
    }

    if (! next_p)
      next_p = &post_last;

    int next_color_pos = (int)(next_p->pos * n_colors) + 1;

    if (next_color_pos <= color_pos)
      next_color_pos = color_pos + 1;

    for (; color_pos < next_color_pos; color_pos ++) {
      float prevpos = p.pos;
      float nextpos = next_p->pos;
      float currentpos = ((float)color_pos) / n_colors;
      float fade;
      if ((nextpos - prevpos) < 1e-3)
        fade = 0.5;
      else
        fade = (currentpos - prevpos) / (nextpos - prevpos);
      float rfade = 1.0 - fade;
      float r = rfade * p.r  +  fade * next_p->r;
      float g = rfade * p.g  +  fade * next_p->g;
      float b = rfade * p.b  +  fade * next_p->b;

      set_color(palette, color_pos, r, g, b, format);
    }

    p = *next_p;
  }
}

void render(SDL_Surface *screen, const int winW, const int winH,
            palette_t *palette, double *pixbuf, const int W, const int H)
{   
  // Lock surface if needed
  if (SDL_MUSTLOCK(screen)) 
    if (SDL_LockSurface(screen) < 0) 
      return;


  int x, y;
  int pitch = screen->pitch / sizeof(Uint32) - winW;
#if 0
  double max = 1e-5;
  double min = 0;
  y = W * H;
  for (x = 0; x < y; x++) {
    double v = pixbuf[x];
    if (v > max)
      max = v;
    if (v < min)
      min = v;
  }
  printf("%f .. %f\n", min, max);
#endif

  Uint32 *screenpos = (Uint32*)(screen->pixels);
  double *pixbufpos;
  for (y = 0; y < H; y++) {
    pixbufpos = pixbuf + y;
    for (x = 0; x < W; x++) {
      int col = *pixbufpos;
      if ((col >= palette->len) || (col < 0)) {
        col = 0;
        *pixbufpos = 0;
      }
      Uint32 raw = palette->colors[(int)col];
      *screenpos = raw;
      screenpos ++;
      pixbufpos += H;
    }
    screenpos += pitch;
  }

#if 0
  {
    int i, l;
    l = palette->len;
    if (l > W*H) {
      l = W*H;
    }
    for (i = 0; i < l; i++) {
      ((Uint32*)screen->pixels)[i] = palette->colors[i];
    }
  }
#endif


  // Unlock if needed
  if (SDL_MUSTLOCK(screen)) 
    SDL_UnlockSurface(screen);

  // Tell SDL to update the whole screen
  SDL_UpdateRect(screen, 0, 0, winW, winH);
}

int main()
{
  // apply FFT to real 2D data.

  double *in;
  double *apex;
  int W = 640;
  int H = 480;
  int half_H = (H / 2) + 1;
  int x;
  int y;
  fftw_complex *out;
  fftw_complex *apex_f;
  fftw_plan plan_backward;
  fftw_plan plan_apex;
  fftw_plan plan_forward;
  unsigned int seed = 123456789;
  srand(seed);

  in = (double *) malloc_check(sizeof(double) * W * H);
  for(x = 0; x < W; x++)
  {
    for(y = 0; y < H; y++)
    {
#if 1
      in[x*H+y] =  ( double ) rand ( ) / ( RAND_MAX );
#else
      in[x*H+y] =  0;
#endif
    }
  }
  in[(H/2) + (W/2)*H] = 1;
  in[(H/2)+3 + (W/2 + 3)*H] = 1;
  in[10 + (20)*H] = 1;
  in[H-3 + (W-3)*H] = 1;

  y = W * H;
  for (x = 0; x < y; x++) {
    in[x] *= PALETTE_LEN -10;
  }


  apex = (double*)malloc_check(sizeof(double) * W * H);
  double apex_sum = 0;
  for(x = 0; x < W; x++)
  {
    for(y = 0; y < H; y++)
    {
      double dist = 0;
      int xx = x;
      int yy = y;
      if (xx >= W/2)
        xx = W - x;
      if (yy >= H/2)
        yy = H - y;
      dist = sqrt(xx*xx + yy*yy);
      double v = 8.01 - dist;
      if (v < 0)
        v = 0;
#if 1
      if (x == 2 && y == 1)
        v = 302.1;
#endif

#if 0
      if (x == W / 2 && y == H / 2)
        v = 850;
#endif

#if 0
      if (x < W/2 || y > H / 2)
        v = -v * 1.85;
#endif
#if 0
      if (x == W/3-1 && y == H/3-1)
        v = 200;
      if (x == W/3 && y == H/3)
        v = -200;
#endif
      apex_sum += v;
      apex[x*H+y] = v;
    }
  }

  double burn = 1.005;
  double apex_mul = (burn / (W*H)) / apex_sum;
  printf("%f %f\n", apex_sum, apex_mul);

  y = W * H;
  for (x = 0; x < y; x++) {
    apex[x] *= apex_mul;
  }

  apex_f = fftw_malloc(sizeof(fftw_complex) * W * half_H);
  plan_apex = fftw_plan_dft_r2c_2d(W, H, apex, apex_f, FFTW_ESTIMATE);
  fftw_execute(plan_apex);


  out = fftw_malloc(sizeof(fftw_complex) * W * half_H);
  plan_forward = fftw_plan_dft_r2c_2d(W, H, in, out, FFTW_ESTIMATE);
  plan_backward = fftw_plan_dft_c2r_2d(W, H, out, in, FFTW_ESTIMATE);



  SDL_Surface *screen;

  if ( SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) < 0 ) 
  {
    fprintf(stderr, "Unable to init SDL: %s\n", SDL_GetError());
    exit(1);
  }
    
  int winW = W;
  int winH = H;
  screen = SDL_SetVideoMode(winW, winH, 32, SDL_SWSURFACE);
  if ( screen == NULL ) 
  {
    fprintf(stderr, "Unable to set %dx%d video: %s\n", winW, winH, SDL_GetError());
    exit(1);
  }
  
  SDL_WM_SetCaption("burnscope", "burnscope");
  SDL_ShowCursor(SDL_DISABLE);

#if 0
#define n_palette_points 2
  palette_point_t palette_points[n_palette_points] = {
    { 0., 0, 0, 0 },
    { 1., 1, 1, 1 },
  };
#else 
#define n_palette_points 11
  palette_point_t palette_points[n_palette_points] = {
    { 0./6, 1, 1, 1 },
    { 0.5/6, 1, .9, 0 },
    { 1./6, 1, .1, 1 },
    { 1.5/6, 0, 0, 1 },
    { 3./6, .5, 0, .7 },
    { 3.5/6, 0, 1, .7 },
    { 4.5/6, .2, .8, .2 },
    { 4.8/6, 0, 0, 1 },
    { 5.25/6, .8, .8, 0 },
    { 5.55/6, .8, .2, 0.4 },
    { 5.85/6, .0,.60,.50 },
  };
#endif

  palette_t palette;
  make_palette(&palette, PALETTE_LEN,
               palette_points, n_palette_points,
               screen->format);

  bool running = true;
  int frame_period = 50;
  int last_ticks = SDL_GetTicks() - frame_period;

  while (running)
  {
    bool do_render = false;

    int elapsed = SDL_GetTicks() - last_ticks;
    if (elapsed > frame_period) {
      last_ticks += frame_period * (elapsed / frame_period);
      do_render = true;
    }

    if (do_render) {
      render(screen, winW, winH, &palette, in, W, H);
      fftw_execute(plan_forward);

#if 1
      for (x = 0; x < W; x++) {
        for (y = 0; y < half_H; y++) {
          double *o = out[x*half_H + y];
          double *af = apex_f[x*half_H + y];
          double a, b, c, d;
          a = o[0]; b = o[1];
          c = af[0]; d = af[1];
#if 1
          o[0] = (a*c - b*d);
          o[1] = (b*c + a*d);
#else
          double l = sqrt(c*c + d*d);
          o[0] *= l;
          o[1] *= l;
#endif
        }
      }
#endif

      fftw_execute(plan_backward);
    }
    else
      SDL_Delay(5);

    SDL_Event event;
    while (SDL_PollEvent(&event)) 
    {
      switch (event.type) 
      {
        case SDL_KEYDOWN:
          // If escape is pressed, return (and thus, quit)

          switch(event.key.keysym.sym) {
            case SDLK_ESCAPE:
              running = false;
              break;

            default:
              break;
          }
          break;

        case SDL_QUIT:
          running = false;
          break;
      }
    }

  }
  SDL_Quit();

  fftw_destroy_plan(plan_apex);
  fftw_destroy_plan(plan_forward);
  fftw_destroy_plan(plan_backward);

  free(in);
  free(apex);
  fftw_free(out);
  fftw_free(apex_f);

  return 0;
}
