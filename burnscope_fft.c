/* burnscope.c
 * (c) 2014 Neels Hofmeyr <neels@hofmeyr.de>
 *
 * This file is part of burnscope, published under the GNU General Public
 * License v3.
 */

#include <fftw3.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <SDL/SDL.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <stdint.h>

#include "dope.h"
#include "shit.h"
#include "human.h"
#include "kill.h"

#define PALETTE_LEN_BITS 12
#define PALETTE_LEN (1 << PALETTE_LEN_BITS)
#define SEED_VAL (0.5 * PALETTE_LEN)

typedef double pixel_t;

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

typedef enum {
  symm_none = 0,
  symm_x = 1,
  symm_y = 2,
  symm_xy = 3,
  symm_point = 4
} symmetry_t;
#define SYMMETRY_KINDS 5


char *symmetry_name[SYMMETRY_KINDS] = {
    "asymmetrical",
    "x-symmetrical (about vertical axis)",
    "y-symmetrical (about horizontal axis)",
    "x- and y-symmetrical (about vertical and horizontal axes)",
    "point-symmetrical"
  };

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


#define min(A,B) ((A) > (B)? (B) : (A))
#define max(A,B) ((A) > (B)? (A) : (B))


int W = 320;
int H = 200;
pixel_t *pixbuf = NULL;
pixel_t *apex = NULL;
fftw_complex *pixbuf_f;
fftw_plan plan_backward;
fftw_plan plan_forward;
fftw_complex *apex_f;
fftw_plan plan_apex;

void make_apex(double apex_r, double burn_factor);

void fft_init(void) {
  int x;
  int y;

  pixbuf = (double *) malloc_check(sizeof(double) * W * H);
  for(x = 0; x < W*H; x++) {
#if 1
      pixbuf[x] =  ( double ) rand ( ) / ( RAND_MAX );
#else
      pixbuf[x] =  0;
#endif
  }
  pixbuf[(H/2) + (W/2)*H] = 1;
  pixbuf[(H/2)+3 + (W/2 + 3)*H] = 1;
  pixbuf[10 + (20)*H] = 1;
  pixbuf[H-3 + (W-3)*H] = 1;

  y = W * H;
  for (x = 0; x < y; x++) {
    pixbuf[x] *= PALETTE_LEN -10;
  }

  int half_W = (W / 2) + 1;
  apex = (double*)malloc_check(sizeof(double) * H * W);
  apex_f = fftw_malloc(sizeof(fftw_complex) * H * half_W);
  plan_apex = fftw_plan_dft_r2c_2d(H, W, apex, apex_f, FFTW_ESTIMATE);

  pixbuf_f = fftw_malloc(sizeof(fftw_complex) * H * half_W);
  plan_forward = fftw_plan_dft_r2c_2d(H, W, pixbuf, pixbuf_f, FFTW_ESTIMATE);
  plan_backward = fftw_plan_dft_c2r_2d(H, W, pixbuf_f, pixbuf, FFTW_ESTIMATE);

  make_apex(8.01, 1.005);
}

void fft_destroy(void) {
  fftw_destroy_plan(plan_apex);
  fftw_destroy_plan(plan_forward);
  fftw_destroy_plan(plan_backward);
  free(pixbuf);
  free(apex);
  fftw_free(pixbuf_f);
  fftw_free(apex_f);
  pixbuf = NULL;
  apex = NULL;
}

void make_apex(double apex_r, double burn_factor) {
  int x, y;

  double apex_sum = 0;
  double apex_r2 = apex_r * apex_r;
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
      dist = xx*xx + yy*yy;
      double v = apex_r2 - dist;
      if (v < 0)
        v = 0;
#if 0
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
      apex[x+y*W] = v;
    }
  }

  double apex_mul = (burn_factor / (W*H)) / apex_sum;

  y = W * H;
  for (x = 0; x < y; x++) {
    apex[x] *= apex_mul;
  }
  fftw_execute(plan_apex);
}


void burn(void) {
  int x;
  int y;
  int half_W = (W / 2) + 1;
  fftw_execute(plan_forward);

  // complex multiplication --> convolution of pixbuf with apex.
  for (x = 0; x < H*half_W; x++) {
    double *pf = pixbuf_f[x];
    double *af = apex_f[x];
    double a, b, c, d;
    a = pf[0]; b = pf[1];
    c = af[0]; d = af[1];
    pf[0] = (a*c - b*d);
    pf[1] = (b*c + a*d);
  }

  fftw_execute(plan_backward);
}


void mirror_x(pixel_t *pixbuf, const int W, const int H) {
  int x, y;
  int x_fold = W >> 1;
  pixel_t *pos_to, *pos_from;

  pos_from = pixbuf + x_fold - 1;
  pos_to = pixbuf + W - x_fold;
  int pitch_to = W - x_fold;
  int pitch_from = W + x_fold;

  for (y = 0; y < H; y ++) {
    for (x = x_fold; x < W; x ++) {
      *(pos_to++) = *(pos_from--);
    }
    pos_from += pitch_from;
    pos_to += pitch_to;
  }
}

void mirror_y(pixel_t *pixbuf, const int W, const int H) {
  int x;
  int y_fold = H >> 1;
  pixel_t *pos_to, *pos_from, *end;
  end = pixbuf + W * H;

  pos_from = pixbuf + (y_fold-1) * W;
  pos_to = pixbuf + (H - y_fold) * W;

  int pitch_from = -2 * W;

  while (pos_to < end) {
    for (x = 0; x < W; x++)
      *(pos_to++) = *(pos_from++);
    pos_from += pitch_from;
  }
}

void mirror_p(pixel_t *pixbuf, const int W, const int H) {
  int x;
  int y_fold = (H >> 1) + (H & 1);
  pixel_t *pos_to, *pos_from, *end;
  end = pixbuf + W * H;

  pos_from = pixbuf + (y_fold-1) * W + (W - 1);
  pos_to = pixbuf + (H - y_fold) * W;

  while (pos_to < end) {
    for (x = 0; x < W; x++)
      *(pos_to++) = *(pos_from--);
  }
}


void render(SDL_Surface *screen, const int winW, const int winH,
            palette_t *palette, pixel_t *pixbuf, 
            int multiply_pixels, int colorshift)
{   
  // Lock surface if needed
  if (SDL_MUSTLOCK(screen)) 
    if (SDL_LockSurface(screen) < 0) 
      return;

  assert((W * multiply_pixels) == winW);
  assert((H * multiply_pixels) == winH);

  int x, y;
  int mx, my;
  int pitch = screen->pitch / sizeof(Uint32) - winW;


  Uint32 *screenpos = (Uint32*)(screen->pixels);
  pixel_t *pixbufpos = pixbuf;
  int winy = 0;
  for (y = 0; y < H; y++) {
    pixel_t *pixbuf_y_pos = pixbufpos;

    for (my = 0; my < multiply_pixels; my ++, winy++) {
      pixbufpos = pixbuf_y_pos;
      int winx = 0;
      for (x = 0; x < W; x++) {
        pixel_t pix = *pixbufpos;
        if (pix >= palette->len) {
          pix -= palette->len;
          *pixbufpos = pix;
        }
        unsigned int col = (unsigned int)pix + colorshift;
        col %= palette->len;
        
        Uint32 raw = palette->colors[col];
        
        for (mx = 0; mx < multiply_pixels; mx++, winx++) {
          *screenpos = raw;
          screenpos ++;
        }
        pixbufpos ++;
      }
      screenpos += pitch;
    }
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

void seed1(pixel_t *pixbuf, const int W, const int H, int x, int y,
           pixel_t val) {
  if ((x < 0) || (x >= W) || (y < 0) || (y >= H))
    return;
  pixbuf[x + y * W] += val;
}

void seed(pixel_t *pixbuf, const int W, const int H, int x, int y,
          pixel_t val, int apex_r) {
  int rx, ry;
  for (ry = -apex_r; ry <= apex_r; ry++) {
    for (rx = -apex_r; rx <= apex_r; rx++) {
      seed1(pixbuf, W, H, x + rx, y + ry, val);
    }
  }
}

void seed_image(int x, int y, char *img, int w, int h) {
  pixel_t *pixbuf_pos = pixbuf + y * W + x;
  pixel_t *pixbuf_end = pixbuf + W * H;
  int pixbuf_pitch = W - w;

  char *img_pos = img;
  int xx, yy;
  for (yy = 0; yy < h; yy++) {
    for (xx = 0; (xx < w) && (pixbuf_pos < pixbuf_end); xx++) {
      if (*img_pos) {
        *pixbuf_pos = SEED_VAL;
      }
      img_pos ++;
      pixbuf_pos ++;
    }
    pixbuf_pos += pixbuf_pitch;
  }
}


int main(int argc, char *argv[])
{
  int multiply_pixels = 2;
  double apex_r = 8.01;
  double burn_factor = 1.005;
  int frame_period = 70;
  bool usage = false;
  bool error = false;
  bool start_blank = false;
  symmetry_t symm = symm_none;

  int c;
  int random_seed = time(NULL);

  char *out_stream_path = NULL;
  FILE *out_stream = NULL;

  while (1) {
    c = getopt(argc, argv, "a:g:m:p:r:u:O:P:ABh");
    if (c == -1)
      break;
   
    switch (c) {
      case 'g':
        {
          char arg[strlen(optarg) + 1];
          strcpy(arg, optarg);
          char *ch = arg;
          while ((*ch) && ((*ch) != 'x')) ch ++;
          if ((*ch) == 'x') {
            *ch = 0;
            ch ++;
            W = atoi(arg);
            H = atoi(ch);

          }
          else {
            fprintf(stderr, "Invalid -g argument: '%s'\n", optarg);
            exit(-1);
          }
        }
        break;

      case 'm':
        multiply_pixels = atoi(optarg);
        break;

      case 'p':
        frame_period = atoi(optarg);
        break;

      case 'a':
        apex_r = atof(optarg);
        break;

      case 'u':
        burn_factor = atof(optarg);
        break;

      case 'r':
        random_seed = atoi(optarg);
        break;

      case 'O':
        out_stream_path = optarg;
        break;

      case 'B':
        start_blank = true;
        break;

      case 'A':
        symm = symm_none;
        break;

      case '?':
        error = true;
      case 'h':
        usage = true;
        break;

    }
  }

  if (usage) {
    if (error)
      printf("\n");
    printf(
"burnscope v0.1\n"
"(c) 2014 Neels Hofmeyr <neels@hofmeyr.de>\n"
"Published under the GNU General Public License v3.\n\n"
"Burnscope produces a mesmerizing animation that I discovered by accident when I\n"
"was a teenager. I've recreated it in memories of old times. It repeatedly\n"
"applies a simple underdamped blur algorithm to a seed image, allowing the color\n"
"values to wrap when overflowing. If you can explain how this staggering\n"
"everchanging complexity can spring from such a simple algorithm and just one\n"
"pixel as seed, please send me an email ;)\n"
"\n"
"Usage example:\n"
"  burnscope -g 320x200 -m 2 -p 70\n"
"\n"
"Options:\n"
"\n"
"  -g WxH   Set animation width and height in number of pixels.\n"
"  -p ms    Set frame period to <ms> milliseconds (slow things down).\n"
"           If zero, run as fast as possible. Default is %d.\n"
"  -m N     Multiply each pixel N times in width and height, to give a larger\n"
"           picture. This will also multiply the window size.\n"
"  -a W     Set apex radius, i.e. the blur distance. Default is %.3f.\n"
"  -u N.n   Set underdampening factor (decimal). Default is %.3f.\n"
"           Reduces normal blur dampening by this factor.\n"
"  -r seed  Supply a random seed to start off with.\n"
"  -B       Start out blank. (Use 's' key to plant seeds while running.)\n"
, frame_period, apex_r, burn_factor
);
    if (error)
      return 1;
    return 0;
  }

  const int maxpixels = 1e4;

  if ((W < 3) || (W > maxpixels) || (H < 3) || (H > maxpixels)) {
    fprintf(stderr, "width and/or height out of bounds: %dx%d\n", W, H);
    exit(-1);
  }

  {
    double was_apex_r = apex_r;
    int max_dim = max(W, H);
    apex_r = min(max_dim, apex_r);
    apex_r = max(1, apex_r);
    if (apex_r != was_apex_r) {
      fprintf(stderr, "Invalid apex radius (-a %f). Forcing %f.", was_apex_r, apex_r);
    }
  }

  float minuscule = 1e-3;
  if ((burn_factor > -minuscule) && (burn_factor < minuscule)) {
    fprintf(stderr, "Underdampening too close to zero (-u). Limit is %f.\n",
        minuscule);
    exit(-1);
  }

  int winW = W;
  int winH = H;

  if (multiply_pixels > 1) {
    winW *= multiply_pixels;
    winH *= multiply_pixels;
  }
  else
    multiply_pixels = 1;

  if ( (winW > maxpixels) || (winH > maxpixels) ) {
    fprintf(stderr, "pixel multiplication is too large: %dx%d times %d = %dx%d\n",
            W, H, multiply_pixels, winW, winH);
    exit(-1);
  }

  if (out_stream_path) {
    out_stream = fopen(out_stream_path, "w");
#if 0
    Uint32 ww = winW;
    Uint32 hh = winH;
    fwrite(&ww, sizeof(ww), 1, out_stream);
    fwrite(&hh, sizeof(hh), 1, out_stream);
#endif
  }


  SDL_Surface *screen;

  if ( SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) < 0 ) 
  {
    fprintf(stderr, "Unable to init SDL: %s\n", SDL_GetError());
    exit(1);
  }
    
  screen = SDL_SetVideoMode(winW, winH, 32, SDL_SWSURFACE);
  if ( screen == NULL ) 
  {
    fprintf(stderr, "Unable to set %dx%d video: %s\n", winW, winH, SDL_GetError());
    exit(1);
  }
  
  SDL_WM_SetCaption("burnscope", "burnscope");
  SDL_ShowCursor(SDL_DISABLE);

#if 1
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
#else
#define n_palette_points 2
  palette_point_t palette_points[n_palette_points] = {
    { 0, 0, 0, 0 },
//    { 0.5, 0,0,0 },
    { 0.5 + 3./256, 0, .8, 0 },
  //  { 0.5 + 6./256, 0, .0, 0 },
  };
#endif

  palette_t palette;
  make_palette(&palette, PALETTE_LEN,
               palette_points, n_palette_points,
               screen->format);

  printf("random seed: %d\n", random_seed);
  srandom(random_seed);

  fft_init();

  if (! start_blank) {
    int i, j;
    j = 2*apex_r + 1;
    j *= j;
    j = W * H / j;
    for (i = 0; i < j; i ++) {
      seed(pixbuf, W, H, random() % (W), random() % (H), SEED_VAL, apex_r);
    }
  }

  if (0)
  {
    int x;
    for (x = 0; x < (W*H); x++)
      printf("%f ", pixbuf[x]);
  }

  int last_ticks = SDL_GetTicks() - frame_period;
  bool seed_key_down = false;
  int do_seed = 0;

  int frames_rendered = 0;
  float wavy = 0;
  float wavy_amp = .006;
  int colorshift = 0;
  bool running = true;
  bool do_stop = false;
  bool do_go = false;
  bool do_wavy = false;
  bool do_stutter = false;

  while (running)
  {
    bool do_render = false;

    float t = (float)frames_rendered / 100.;
    wavy = sin(t);
    colorshift = palette.len * (0.5 + 0.5 * cos(t*M_PI/50));
    
    float use_burn = burn_factor;
    if (do_wavy)
      use_burn += (wavy_amp * wavy);

    if (frame_period < 1)
      do_render = true;
    else {
      int elapsed = SDL_GetTicks() - last_ticks;
      if (elapsed > frame_period) {
        last_ticks += frame_period * (elapsed / frame_period);
        do_render = true;
      }
    }


    if (do_render) {
      static bool stopped = false;
      if (do_stop) {
        stopped = true;
        do_stop = false;
        // and do just one rendering
      }
      else
      if (do_go) {
        stopped = false;
        do_go = false;
      }
      else
      if (stopped)
        do_render = false;

      bool do_calc = true;

      if (do_stutter) {
        static char stutter_count = 0;
        if ((stutter_count++) > 2)
          stutter_count = 0;
        else
          do_calc = false;
      }

      if (do_calc) {
        if (seed_key_down) {
          static int seed_slew = 0;
          if ((++ seed_slew) > 1) {
            seed_slew = 0;
            do_seed ++;
          }
        }

        while (do_seed) {
          do_seed --;
          int seedx = random() % W;
          int seedy = random() % H;
          seed(pixbuf, W, H, seedx, seedy, SEED_VAL, apex_r);
          if ((symm == symm_x) || (symm == symm_xy))
            seed(pixbuf, W, H, W - seedx, seedy, SEED_VAL, apex_r);
          if ((symm == symm_y) || (symm == symm_xy))
            seed(pixbuf, W, H, seedx, H - seedy, SEED_VAL, apex_r);
          if (symm == symm_point)
            seed(pixbuf, W, H, W - seedx, H - seedy, SEED_VAL, apex_r);
        }


        burn();

        if (symm == symm_x)
          mirror_x(pixbuf, W, H);
        else
        if (symm == symm_xy)
          mirror_x(pixbuf, W, H);
        if ((symm == symm_y) || (symm == symm_xy))
          mirror_y(pixbuf, W, H);
        if (symm == symm_point)
          mirror_p(pixbuf, W, H);
      }

      if (do_render) {
        render(screen, winW, winH, &palette, pixbuf, multiply_pixels, colorshift);
      }

      if (out_stream) {
        fwrite(screen->pixels, sizeof(Uint32), winW * winH, out_stream);
      }

      frames_rendered ++;

      static double was_apex_r = 0;
      static double was_burn = 0;

      if ((was_apex_r != apex_r) || (was_burn != use_burn)) {
        make_apex(apex_r, use_burn);
        was_apex_r = apex_r;
        was_burn = use_burn;
      }

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


          {
          int c = event.key.keysym.sym;
          switch(c) {
            case SDLK_ESCAPE:
              running = false;
              break;

            case SDLK_RIGHT:
              burn_factor += .0002;
              break;
            case SDLK_LEFT:
              burn_factor -= .0002;
              break;
            case SDLK_UP:
              wavy_amp += .0001;
              break;
            case SDLK_DOWN:
              wavy_amp -= .0001;
              break;

            case 's':
              do_seed ++;
              seed_key_down = true;
              break;

            case 'b':
              bzero(pixbuf, W * H * sizeof(pixel_t));
              break;

            case 'm':
              symm = (symm + 1) % SYMMETRY_KINDS;
              break;

            case 'q':
              burn_factor -= .002;
              break;
            case 'w':
              burn_factor -= .0003;
              break;
            case 'e':
              burn_factor = (burn_factor + burn_factor + 1.005) / 3.;
              break;
            case 'r':
              burn_factor += .0003;
              break;
            case 't':
              burn_factor += .002;
              break;

            case '`':
              do_wavy = ! do_wavy;
              break;

            case '-':
              apex_r = max(0.5, apex_r / 1.1);
              break;

            case '+':
            case '=':
              apex_r = min(W, apex_r * 1.1);
              break;

            case '/':
              do_go = true;
              do_stutter = false;
              break;

            case '.':
              do_stop = true;
              break;

            case ',':
              do_stutter = ! do_stutter;
              do_go = true;
              break;

#define drop_img(name)  seed_image(random() % (30 + W - name##_width), random() %(30 + H- name##_height), name##_data, name##_width, name##_height)

            case 'u':
              drop_img(human);
              break;

            case 'i':
              drop_img(kill);
              break;

            case 'p':
              drop_img(dope);
              break;

            case 'o':
              drop_img(shit);
              break;


            default:

              if ((c >= '1') && (c <= '9')) {
                apex_r = 1 + c - '1';
              }
              break;
          }
          }

          printf("burn=%f  wavy=%s wavy_amp=%f  symm=%s  apex_r=%f  stutter=%s\n",
            burn_factor, do_wavy? "on" : "off", wavy_amp, symmetry_name[symm], apex_r, do_stutter? "on":"off");
          break;

        case SDL_KEYUP:
          if (event.key.keysym.sym == 's') {
            seed_key_down = false;
          }
          break;

        case SDL_QUIT:
          running = false;
          break;
      }
    }

  }

  printf("\n");
  printf("%d frames rendered\n", frames_rendered);
  if (out_stream) {
    fclose(out_stream);
    out_stream = NULL;
    printf("suggestion:\n"
        "ffmpeg -vcodec rawvideo -f rawvideo -pix_fmt rgb32 -s %dx%d -i %s  -vcodec libx264 -b 20000k %s.avi\n", W * multiply_pixels, H * multiply_pixels, out_stream_path, out_stream_path);
  }
  fft_destroy();
  SDL_Quit();
  return 0;
}
