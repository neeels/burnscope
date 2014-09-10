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
#include <SDL/SDL_thread.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <stdint.h>

#define PALETTE_LEN_BITS 12
#define PALETTE_LEN (1 << PALETTE_LEN_BITS)
#define SEED_VAL (0.5 * PALETTE_LEN)

#define min(A,B) ((A) > (B)? (B) : (A))
#define max(A,B) ((A) > (B)? (A) : (B))


typedef double pixel_t;

#include "images.h"

int n_images = 0;
image_t *images = NULL;

const float minuscule = 1e-3;

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



int W = 1920 / 3;
int H = 1080 / 3;
int min_W_H, max_W_H;
pixel_t *pixbuf = NULL;
pixel_t *apex = NULL;
fftw_complex *pixbuf_f;
fftw_plan plan_backward;
fftw_plan plan_forward;
fftw_complex *apex_f;
fftw_plan plan_apex;

void make_apex(double apex_r, double burn_factor, char apex_opt);

void fft_init(void) {
  int x;
  int y;

  fftw_init_threads();
  fftw_plan_with_nthreads(2);

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

  bzero(apex, W*H*sizeof(pixel_t));
  make_apex(8.01, 1.005, 0);
}

void fft_destroy(void) {
  fftw_cleanup_threads();
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

void make_apex(double apex_r, double burn_factor, char apex_opt) {
  int x, y;

  apex_r = min(apex_r, min_W_H/2 - 2);

  if (fabs(burn_factor) < minuscule) {
    if (burn_factor < 0)
      burn_factor = -minuscule;
    else
      burn_factor = minuscule;
  }
  
  double apex_sum = 0;
  double apex_r2 = apex_r * apex_r;
  int apex_r_i = apex_r;

  static int last_apex_r = 0;
  int overwrite_r = max(apex_r_i, last_apex_r);

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

      double v;
      if ((xx > apex_r_i) || (yy > apex_r_i))
        v = 0;
      else
      {
        dist = xx*xx + yy*yy;
        v = apex_r2 - dist;
        if (v < 0)
          v = 0;
      }

#if 1
      if (apex_opt)
      switch(apex_opt) {
        default:
          break;

        case 1:
          if (x > W/2 || y < H / 2)
            v = -v * 1.85;
          break;

        case 2:
          if (x < W/2 || y < H / 2)
            v = -v * 1.85;
          break;

        case 3:
          if (x < W/2 || y > H / 2)
            v = -v * 1.85;
          break;

        case 4:
          if (x > W/2 || y > H / 2)
            v = -v * 1.85;
          break;
      }
#endif


      apex_sum += v;
      apex[x+y*W] = v;

      if (y == overwrite_r)
        y = H - overwrite_r - 1;
    }
    if (x == overwrite_r)
      x = W - overwrite_r - 1;
  }
#if 0
  for(x = 0; x < apex_r; x++)
  {
    for(y = 0; y < apex_r; y++)
    {
      v = apex_r2 - (x*x + y*y);
      if (v < 0)
        v = 0;
      apex[x+y*W] = v;
      apex_sum += v;
      if ((x > 0) && (y > 0)) {
        apex[(W - x)+y*W] = v;
        apex_sum += v;
        apex[(W - x)+(H - y)*W] = v;
        apex_sum += v;
        apex[x+(H - y)*W] = v;
        apex_sum += v;
      }

    }
  }
#endif

#if 0
  if (apex_opt) {
    int at = 0;
    switch(apex_opt) {
      default:
        break;

      case 1:
        at = W / 3 + W*(H/3);
        break;

      case 2:
        at = W - apex_r;
        break;

      case 3:
        at = (H - apex_r)*W;
        break;

      case 4:
        at = (H- apex_r)*W + W - apex_r;
        break;

    }

    if (at) {
      for (y = 0; y < apex_r; y ++, at+=(W - apex_r))
      {
        for (x = 0; x < apex_r; x ++, at++) {
          double was = apex[at];
          apex[at] *= -1.85;
          apex_sum += apex[at] - was;
        }
      }
    }
  }
#endif

  double apex_mul = (burn_factor / (W*H)) / apex_sum;

  y = W * H;
  for (x = 0; x < y; x++) {
    apex[x] *= apex_mul;
  }
  fftw_execute(plan_apex);
  last_apex_r = apex_r_i;
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
      pixel_t v = min(*pos_to, *pos_from);
      *pos_to = v;
      *pos_from = v;
      pos_to++;
      pos_from--;
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
    for (x = 0; x < W; x++) {
      pixel_t v = min(*pos_to, *pos_from);
      *pos_to = v;
      *pos_from = v;
      pos_to++;
      pos_from++;
    }
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
    for (x = 0; x < W; x++) {
      pixel_t v = min(*pos_to, *pos_from);
      *pos_to = v;
      *pos_from = v;
      pos_to++;
      pos_from--;
    }
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
          pix -= palette->len * (int)(pix) / palette->len;
          *pixbufpos = pix;
        }
        else
        if (pix < 0.001) {
          pix += palette->len * (1 + (int)(-pix) / palette->len);
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

void seed_image(int x, int y, pixel_t *img, int w, int h) {
  pixel_t *pixbuf_pos = pixbuf + y * W + x;
  pixel_t *pixbuf_end = pixbuf + W * H;
  int pixbuf_pitch = max(0, W - w);

  pixel_t *img_pos = img;
  int xx, yy;
  for (yy = 0; yy < h; yy++) {
    for (xx = 0; (xx < w) && (pixbuf_pos < pixbuf_end); xx++) {
      (*pixbuf_pos) += (*img_pos) * (PALETTE_LEN >> 1);
      img_pos ++;
      pixbuf_pos ++;
    }
    pixbuf_pos += pixbuf_pitch;
  }
}

volatile bool running = true;
volatile int frames_rendered = 0;

volatile int avg_frame_period = 0;
#define AVG_SHIFTING 3
int want_frame_period = 40;

SDL_Surface *screen;
int winW;
int winH;
palette_t palette;
FILE *out_stream = NULL;
FILE *out_params = NULL;
FILE *in_params = NULL;
int multiply_pixels = 2;
int colorshift;

typedef struct {
  int random_seed;
} init_params_t;

typedef struct {
  float apex_r;
  char apex_opt;
  float burn_factor;
  float axis_colorshift;
  float axis_seed;
  bool do_stop;
  bool do_go;
  bool do_wavy;
  bool do_stutter;
  bool do_blank;
  bool do_maximize;
  bool force_symm;
  symmetry_t symm;
  float seed_r;
  int n_seed;
  float wavy_amp;
  int please_drop_img;
} params_t;

const int params_file_id = 0x23315;
const int params_version = 1;

init_params_t ip;
params_t p = { 8.01, 0, 1.002, 0., 0., false, false, false, false, false, false, false, symm_none, 3, 0, .006, -1};


int normalize_colorshift = 0;

void maximize(void) {
  pixel_t *pos;
  pixel_t *end = pixbuf + (W * H);
  pixel_t max_val = -1;
  for (pos = pixbuf; pos < end; pos++) {
    max_val = max(max_val, *pos);
  }
  pixel_t diff = (pixel_t)PALETTE_LEN - max_val;

  for (pos = pixbuf; pos < end; pos++) {
    *pos += diff;
  }

  normalize_colorshift -= diff;
  printf("normalized %+.2f\n", diff);
}


SDL_sem *please_render;
SDL_sem *please_save;
SDL_sem *rendering_done;
SDL_sem *saving_done;

int render_thread(void *arg) {

  int last_ticks = SDL_GetTicks() - want_frame_period;

  for (;;) {

    SDL_SemWait(please_render);
    if (! running)
      break;

    while (want_frame_period) {
      int elapsed = SDL_GetTicks() - last_ticks;
      if (elapsed >= want_frame_period) {
        last_ticks += want_frame_period * (elapsed / want_frame_period);
        break;
      }
      SDL_Delay(want_frame_period - elapsed);
    }

    if (out_stream) {
      SDL_SemWait(saving_done);
    }

    render(screen, winW, winH, &palette, pixbuf, multiply_pixels, colorshift);

    int t = SDL_GetTicks();

    if (out_stream) {
      SDL_SemPost(please_save);
    }

    frames_rendered ++;
    SDL_SemPost(rendering_done);

    {
      static int last_ticks2 = 0;
      int elapsed = t - last_ticks2;
      last_ticks2 = t;
      avg_frame_period -= avg_frame_period >>AVG_SHIFTING;
      avg_frame_period += elapsed;
    }

  }

  return 0;
}

int save_thread(void *arg) {

  for (;;) {
    SDL_SemWait(please_save);
    if (! running)
      break;

    if (out_stream) {
      fwrite(screen->pixels, sizeof(Uint32), winW * winH, out_stream);
    }
    SDL_SemPost(saving_done);
  }

  return 0;
}

int main(int argc, char *argv[])
{
  bool usage = false;
  bool error = false;
  bool start_blank = false;

  int c;

  char *out_stream_path = NULL;
  char *out_params_path = NULL;
  char *in_params_path = NULL;

  ip.random_seed = time(NULL);

  while (1) {
    c = getopt(argc, argv, "a:g:m:p:r:u:i:o:O:P:ABh");
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
        want_frame_period = atoi(optarg);
        break;

      case 'a':
        p.apex_r = atof(optarg);
        break;

      case 'u':
        p.burn_factor = atof(optarg);
        break;

      case 'r':
        ip.random_seed = atoi(optarg);
        break;

      case 'O':
        out_stream_path = optarg;
        break;

      case 'o':
        out_params_path = optarg;
        break;

      case 'i':
        in_params_path = optarg;
        break;

      case 'B':
        start_blank = true;
        break;

      case 'A':
        p.symm = symm_none;
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
"           Default is '-g %dx%d'.\n"
"  -p ms    Set frame period to <ms> milliseconds (slow things down).\n"
"           If zero, run as fast as possible. Default is %d.\n"
"  -m N     Multiply each pixel N times in width and height, to give a larger\n"
"           picture. This will also multiply the window size.\n"
"  -a W     Set apex radius, i.e. the blur distance. Default is %.3f.\n"
"  -u N.n   Set underdampening factor (decimal). Default is %.3f.\n"
"           Reduces normal blur dampening by this factor.\n"
"  -r seed  Supply a random seed to start off with.\n"
"  -B       Start out blank. (Use 's' key to plant seeds while running.)\n"
"  -O file  Write raw video data to file (grows large quickly). Can be\n"
"           converted to a video file using e.g. ffmpeg.\n"
"  -o file  Write live control parameters to file for later playback, see -i.\n"
"  -i file  Play back previous control parameters (possibly in a different\n"
"           resolution and streaming video to file...)\n"
, W, H, want_frame_period, p.apex_r, p.burn_factor
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

  min_W_H = min(W, H);
  max_W_H = max(W, H);

  {
    double was_apex_r = p.apex_r;
    p.apex_r = min(max_W_H, p.apex_r);
    p.apex_r = max(1, p.apex_r);
    if (p.apex_r != was_apex_r) {
      fprintf(stderr, "Invalid apex radius (-a %f). Forcing %f.", was_apex_r, p.apex_r);
    }
  }

  if ((p.burn_factor > -minuscule) && (p.burn_factor < minuscule)) {
    fprintf(stderr, "Underdampening too close to zero (-u). Limit is %f.\n",
        minuscule);
    exit(-1);
  }

  winW = W;
  winH = H;

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

  printf("%dx%d --> %dx%d\n", W, H, winW, winH);

  if (out_stream_path) {
    out_stream = fopen(out_stream_path, "w");
  }

  if (out_params_path) {
    out_params = fopen(out_params_path, "w");
  }

  if (in_params_path) {
    in_params = fopen(in_params_path, "r");
  }


  int in_params_framelen = sizeof(p);
  int in_params_read_framelen = sizeof(p);

  if (in_params) {
    int v;
    fread(&v, sizeof(v), 1, in_params);
    if (v != params_file_id) {
      fprintf(stderr,
              "This does not appear to be a burnscope parameters file: %s\n",
              in_params_path);
      exit(1);
    }
    fread(&v, sizeof(v), 1, in_params);
    if (v != params_version) {
      printf("Parameter file has different version: %d. I am at %d.\n", v, params_version);
    }
    int ip_len;
    fread(&ip_len, sizeof(ip_len), 1, in_params);
    if (ip_len < sizeof(ip)) {
      printf("Parameter file header is smaller: %d instead of %d\n", ip_len, (int)sizeof(ip));
    }
    if (ip_len > sizeof(ip)) {
      printf("Parameter file header is larger, reading smaller header: %d instead of %d\n", (int)sizeof(ip), ip_len);
      ip_len = sizeof(ip);
    }
    fread(&ip, ip_len, 1, in_params);

    fread(&in_params_framelen, sizeof(in_params_framelen), 1, in_params);
    in_params_framelen = in_params_framelen;
    if (in_params_framelen < sizeof(p)) {
      printf("Parameter file frame size is smaller: %d instead of %d\n", in_params_framelen, (int)sizeof(p));
    }
    if (in_params_framelen > sizeof(p)) {
      printf("Parameter file frame size is larger, reading smaller frame: %d instead of %d\n", (int)sizeof(p), in_params_framelen);
      in_params_read_framelen = sizeof(p);
    }
  }

  printf("random seed: %d\n", ip.random_seed);
  srandom(ip.random_seed);

  if (out_params) {
#define params_write(what) \
    fwrite(&what, sizeof(what), 1, out_params)

    params_write(params_file_id);
    params_write(params_version);
    int l = sizeof(ip);
    params_write(l);
    params_write(ip);
    l = sizeof(p);
    params_write(l);
  }


  read_images("./images", &images, &n_images);


  if ( SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_JOYSTICK) < 0 ) 
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

  const int n_joysticks = SDL_NumJoysticks();

  printf("%d joysticks were found.\n", n_joysticks);

  SDL_Joystick **joysticks = NULL;

  if (n_joysticks) {
    SDL_JoystickEventState(SDL_ENABLE);

    joysticks = malloc(sizeof(SDL_Joystick*) * n_joysticks);
    
    int i;
    for (i = 0; i < n_joysticks; i++)
    {
      printf("%2d: '%s'\n", i, SDL_JoystickName(i));

      SDL_Joystick *j = SDL_JoystickOpen(i);
      printf("    %d buttons  %d axes  %d balls %d hats\n",
             SDL_JoystickNumButtons(j),
             SDL_JoystickNumAxes(j),
             SDL_JoystickNumBalls(j),
             SDL_JoystickNumHats(j)
             );
      joysticks[i] = j;
    }
  }

#if 1
#define n_palette_points 11
  palette_point_t palette_points[n_palette_points] = {
    { 0./6, 1, 1, 1 },
    { 0.5/6, 0.1, 1, .10 },
    { 1./6, 0, .1, .5 },
    { 1.5/6, .3, .3,1 },
    { 3./6, 1, 1, 0 },
    { 3.5/6, 0, 1, 1 },
    { 4.5/6, 0, 0, 1 },
    { 4.8/6, 1, 0, 0 },
    { 5.25/6, 1, 1, 0 },
    { 5.55/6, 1, .6, 0 },
    { 5.85/6, .1,.0,0 },
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

  make_palette(&palette, PALETTE_LEN,
               palette_points, n_palette_points,
               screen->format);

  fft_init();

  if (! start_blank) {
    int i, j;
    j = 2*p.apex_r + 1;
    j *= j;
    j = W * H / j;
    for (i = 0; i < j; i ++) {
      seed(pixbuf, W, H, random() % (W), random() % (H), SEED_VAL, p.apex_r);
    }
  }

  if (0)
  {
    int x;
    for (x = 0; x < (W*H); x++)
      printf("%f ", pixbuf[x]);
  }

  float wavy = 0;
  bool do_print = true;
  float wavy_speed = 3.;
  bool seed_key_down = false;
  double use_burn = 1.002;

  please_render = SDL_CreateSemaphore(0);
  please_save = SDL_CreateSemaphore(0);
  rendering_done = SDL_CreateSemaphore(0);
  saving_done = SDL_CreateSemaphore(1);

  SDL_Thread *render_thread_token = SDL_CreateThread(render_thread, NULL);
  SDL_Thread *save_thread_token = NULL;
  if (out_stream)
    SDL_CreateThread(save_thread, NULL);

  fftw_execute(plan_forward);
  while (running)
  {
    if (in_params) {
      if (! fread(&p, in_params_read_framelen, 1, in_params)) {
        running = false;
      }
      else
      if (in_params_framelen > in_params_read_framelen) {
        // skip trailing bytes if params file frame is larger than my params_t.
        fseek(in_params, in_params_framelen - in_params_read_framelen, SEEK_CUR);
      }
    }

    bool do_calc = true;
    {
      static bool stopped = false;
      if (p.do_stop) {
        stopped = true;
        p.do_stop = false;
        // and do just one rendering
      }
      else
      if (p.do_go) {
        stopped = false;
        p.do_go = false;
      }
      else
      if (stopped)
        do_calc = false;

      if (p.do_stutter) {
        static char stutter_count = 0;
        if ((stutter_count++) > 1)
          stutter_count = 0;
        else
          do_calc = false;
      }

    }

    if (p.do_maximize) {
      //p.do_maximize = false; first save below
      maximize();
    }

    if (p.do_blank) {
      // p.do_blank = false; first save below
      bzero(pixbuf, W * H * sizeof(pixel_t));
    }

    if (do_calc) {
      float t = (float)frames_rendered / 100.;
      wavy = sin(wavy_speed*t);
      
      colorshift = normalize_colorshift;

      {
        static int colorshift_axis_accum = 0;
        colorshift_axis_accum += 70.*p.axis_colorshift*p.axis_colorshift;
        colorshift += colorshift_axis_accum;
      }
      colorshift %= PALETTE_LEN;

      
      use_burn = p.burn_factor;

      if (p.do_wavy) {
        use_burn += (p.wavy_amp * wavy);
      }


      if (seed_key_down) {
        static int seed_slew = 0;
        if ((++ seed_slew) > 1) {
          seed_slew = 0;
          p.n_seed ++;
        }
      }

      if (p.axis_seed) {
        static int seed_slew = 0;
        if ((++ seed_slew) >= 5.*(1. - p.axis_seed)) {
          seed_slew = 0;
          p.n_seed ++;
        }
      }

      // short interruption of if(do_calc) to save parameters frame.
    }

    if (out_params) {
      fwrite(&p, sizeof(p), 1, out_params);
    }
    if (p.do_blank) {
      p.do_blank = false; // blanked above
    }
    if (p.do_maximize) {
      p.do_maximize = false; // maximized above
    }

    if (do_calc) {
      while (p.n_seed) {
        p.n_seed --;
        int seedx = random() % W;
        int seedy = random() % H;
        seed(pixbuf, W, H, seedx, seedy, SEED_VAL, p.seed_r);
        if ((p.symm == symm_x) || (p.symm == symm_xy))
          seed(pixbuf, W, H, W - seedx, seedy, SEED_VAL, p.seed_r);
        if ((p.symm == symm_y) || (p.symm == symm_xy))
          seed(pixbuf, W, H, seedx, H - seedy, SEED_VAL, p.seed_r);
        if (p.symm == symm_point)
          seed(pixbuf, W, H, W - seedx, H - seedy, SEED_VAL, p.seed_r);
      }

      if (p.please_drop_img >= 0 && p.please_drop_img < n_images) {
        image_t *img = &images[p.please_drop_img];
        seed_image(random() % (30 + W - img->width), random() %(30 + H- img->height), img->data, img->width, img->height);
        p.please_drop_img = -1;
      }


      {
        if (p.force_symm) {
          p.force_symm = false;
          if (p.symm == symm_x)
            mirror_x(pixbuf, W, H);
          else
          if (p.symm == symm_xy)
            mirror_x(pixbuf, W, H);
          if ((p.symm == symm_y) || (p.symm == symm_xy))
            mirror_y(pixbuf, W, H);
          if (p.symm == symm_point)
            mirror_p(pixbuf, W, H);
        }
      }

      int x;
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

    SDL_SemPost(please_render);

    {
      static double was_apex_r = 0;
      static double was_burn = 0;
      static char was_apex_opt = 0;

      p.apex_r = fabs(p.apex_r);

      if ((was_apex_r != p.apex_r) || (was_burn != use_burn) || (was_apex_opt != p.apex_opt)) {
        make_apex(p.apex_r, use_burn, p.apex_opt);
        was_apex_r = p.apex_r;
        was_burn = use_burn;
        was_apex_opt = p.apex_opt;
      }

    }

    static bool apex_r_clamp = false;
    static bool burn_clamp = false;

    static char printcount = 0;
    if (printcount++ >= 10) {
      if (printcount >= 50)
        do_print = true;

      if (do_print) {
        do_print = false;
        printcount = 0;
        printf("%.1ffps apex_r=%s%f_opt%d burn=%s%f symm=%s stutter=%s\n",
               1000./(avg_frame_period>>AVG_SHIFTING),
               apex_r_clamp ? "*" : "",
               p.apex_r,p.apex_opt,
               burn_clamp ? "*" : "",
               p.burn_factor,
               symmetry_name[p.symm],
               p.do_stutter? "on":"off");
        fflush(stdout);
      }
    }


    while (running) {

      SDL_Event event;
      while (SDL_PollEvent(&event)) 
      {
        bool update_apex_r = false;
        static double smallaxis_apex_r = 0;
        static double largeaxis_apex_r = 0;
        static double apex_r_center = 10;

        bool update_burn = false;
        static double smallaxis_burn = 0;
        static double largeaxis_burn = 0;
        static double burn_center = 1.001;

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
                p.burn_factor += .0002;
                break;
              case SDLK_LEFT:
                p.burn_factor -= .0002;
                break;
              case SDLK_UP:
                p.wavy_amp += .0001;
                break;
              case SDLK_DOWN:
                p.wavy_amp -= .0001;
                break;

              case ' ':
                p.n_seed ++;
                seed_key_down = true;
                break;

              case 'b':
                p.do_blank = true;
                break;

              case 'm':
                p.symm = (p.symm + 1) % SYMMETRY_KINDS;
                break;

              case '\\':
                p.symm = symm_x;
                p.force_symm = true;
                break;

              case '\'':
                p.symm = symm_point;
                p.force_symm = true;
                break;

              case ';':
                p.symm = symm_none;
                break;

              case 'q':
                p.burn_factor -= .002;
                break;
              case 'w':
                p.burn_factor -= .0003;
                break;
              case 'e':
                p.burn_factor = 1.005;
                p.apex_r = 8.01 * min_W_H / 240.;
                break;
              case 'r':
                p.burn_factor += .0003;
                break;
              case 't':
                p.burn_factor += .002;
                break;

              case '`':
                p.do_wavy = ! p.do_wavy;
                break;

              case '-':
                p.apex_r = max(0.5, p.apex_r / 1.1);
                break;

              case '+':
              case '=':
                p.apex_r = min(W, p.apex_r * 1.1);
                break;

              case '/':
                p.do_go = true;
                p.do_stutter = false;
                break;

              case '.':
                p.do_stop = true;
                break;

              case ',':
                p.do_stutter = ! p.do_stutter;
                p.do_go = true;
                break;

              case 'u':
                p.please_drop_img = 0;
                break;

              case 'i':
                p.please_drop_img = 1;
                break;

              case 'o':
                p.please_drop_img = 2;
                break;

              case 'p':
                p.please_drop_img = 3;
                break;

              case 'a':
                p.apex_opt = 0;
                break;

              case 's':
                p.apex_opt = 1;
                break;

              case 'd':
                p.apex_opt = 2;
                break;

              case 'f':
                p.apex_opt = 3;
                break;

              case 'g':
                p.apex_opt = 4;
                break;

              case '0':
                p.apex_r += (float)min_W_H / 48;
                break;

              case 'l':
                wavy_speed += .5;
                break;

              case 'k':
                wavy_speed -= .5;
                break;

              case '1':
                p.apex_r = 1;
                break;

              default:

                if ((c >= '2') && (c <= '9')) {
                  p.apex_r = ((float)min_W_H / 240.) * (1 + c - '1');
                }
                break;
            }
            }
            do_print = true;
            break;

          case SDL_KEYUP:
            if (event.key.keysym.sym == ' ') {
              seed_key_down = false;
            }
            break;

          case SDL_QUIT:
            running = false;
            break;

          case SDL_JOYAXISMOTION:
            {
              float axis_val = event.jaxis.value;
              axis_val /= 32768;
#if 0
#define axis_clear .2
              if (fabs(axis_val) < axis_clear)
                axis_val = 0;
              else
                axis_val = (axis_val - axis_clear) / (1.0 - axis_clear);
#endif

              //axis_val *= axis_val * axis_val * 1.2;


              switch(event.jaxis.axis) {
                case 3:
                  smallaxis_apex_r = axis_val;
                  update_apex_r = true;
                  do_print = true;
                  break;
                case 4:
                  largeaxis_apex_r = -axis_val;
                  update_apex_r = true;
                  do_print = true;
                  break;
                
                case 0:
                  smallaxis_burn = axis_val;
                  update_burn = true;
                  do_print = true;
                  break;
                case 1:
                  largeaxis_burn = -axis_val;
                  update_burn = true;
                  do_print = true;
                  break;

                case 5:
                  p.axis_seed = (axis_val + .5) / 1.5;
                  if (p.axis_seed < 0)
                    p.axis_seed = 0;
                  p.axis_seed *= p.axis_seed;
                  break;

                case 2:
                  p.axis_colorshift = (1. + axis_val) / 2;
                  break;

                default:
                  printf("axis %d = %.3f\n", event.jaxis.axis, axis_val);
                  break;
              }
            }
            break;

          case SDL_JOYBUTTONDOWN:
            switch (event.jbutton.button) {
              case 4:
                p.do_maximize = true;
                break;

              case 2:
                p.symm = symm_x;
                p.force_symm = true;
                break;

              case 1:
                p.symm = symm_point;
                p.force_symm = true;
                break;

              case 3:
                p.do_blank = true;
                break;

              case 5:
                p.n_seed ++;
                break;

              case 10:
                apex_r_center = p.apex_r;
                apex_r_clamp = true;
                printf("apex_r clamp engage\n");
                break;

              case 9:
                burn_center = p.burn_factor;
                burn_clamp = true;
                break;

              default:
                printf("%2d: button %d = %s\n",
                       event.jbutton.which, event.jbutton.button,
                       event.jbutton.state == SDL_PRESSED ? "pressed" : "released");
                break;
            }
            break;

          case SDL_JOYBUTTONUP:
            switch (event.jbutton.button) {

              default:
                printf("%2d: button %d = %s\n",
                       event.jbutton.which, event.jbutton.button,
                       event.jbutton.state == SDL_PRESSED ? "pressed" : "released");
                break;
            }
            break;

          case SDL_JOYHATMOTION:  /* Handle Hat Motion */
            printf("%2d: hat %d = %d\n",
                   event.jhat.which, event.jhat.hat, event.jhat.value);
            break;

          case SDL_JOYBALLMOTION:  /* Handle Joyball Motion */
            printf("%2d: hat %d += %d, %d\n",
                   event.jball.which, event.jball.ball,
                   event.jball.xrel, event.jball.yrel);
            break;
        }

        #define calc_axis_val(start, center, end, axis_val) \
            ((axis_val <= 0)? \
              ((start) + ((center) - (start)) * (1 + (axis_val))) \
              : \
              ((center) + ((end) - (center)) * (axis_val) * (axis_val)))

        if (update_apex_r) {
          static double largest_val = 0.5;
          double val = (smallaxis_apex_r * .2 + largeaxis_apex_r * .8);
          largest_val = max(largest_val, fabs(val));
          val /= largest_val;
          if (apex_r_clamp && (fabs(val) < .05))
            apex_r_clamp = false;
          if (! apex_r_clamp) {
            double apex_r_min = max(.99, apex_r_center * 0.08);
            double apex_r_max = min(min_W_H/2, apex_r_center * 10);
            p.apex_r = calc_axis_val(apex_r_min, apex_r_center, apex_r_max, val);

            p.seed_r = calc_axis_val(1, 3, min_W_H/5, val);
          }
        }
        if (update_burn) {
          static double largest_val = 0.5;
          double val = (smallaxis_burn * .3 + largeaxis_burn * .7);
          largest_val = max(largest_val, fabs(val));
          val /= largest_val;
          if (burn_clamp && (fabs(val) < .05))
            burn_clamp = false;
          if (! burn_clamp) {
            double burn_min = max(.993, burn_center * (.993/1.0005));
            double burn_max = min(1.025, burn_center * (1.025/1.0005));
            p.burn_factor = calc_axis_val(burn_min, burn_center, burn_max, val);
          }
        }

      }

      if (! running)
        break;

      if (SDL_SemTryWait(rendering_done) == 0)
        break;
      else
        SDL_Delay(5);
    }
  }

  running = false;

  SDL_SemPost(please_render);
  SDL_SemPost(please_render);
  if (out_stream) {
    SDL_SemPost(please_save);
    SDL_SemPost(please_save);
  }
  printf("waiting for render thread...\n");
  SDL_WaitThread(render_thread_token, NULL);
  if (out_stream) {
    printf("waiting for save thread...\n");
    SDL_WaitThread(save_thread_token, NULL);
  }

  SDL_DestroySemaphore(please_render);
  SDL_DestroySemaphore(please_save);
  SDL_DestroySemaphore(rendering_done);
  SDL_DestroySemaphore(saving_done);

  if (joysticks) {
    int i;
    for (i = 0; i < n_joysticks; i++)
    {
      SDL_JoystickClose(joysticks[i]);
    }
    free(joysticks);
    joysticks = NULL;
  }

  printf("\n");
  printf("%d frames rendered\n", frames_rendered);
  if (out_stream) {
    fclose(out_stream);
    out_stream = NULL;
    printf("suggestion:\n"
        "ffmpeg -vcodec rawvideo -f rawvideo -pix_fmt rgb32 -s %dx%d -i %s  -vcodec libx264 -b 20000k %s.avi\n", W * multiply_pixels, H * multiply_pixels, out_stream_path, out_stream_path);
  }
  if (out_params) {
    fclose(out_params);
    out_params = NULL;
  }
  if (in_params) {
    fclose(in_params);
    in_params = NULL;
  }
  fft_destroy();
  SDL_Quit();
  return 0;
}
