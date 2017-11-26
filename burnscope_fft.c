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
#include <SDL2/SDL.h>
#include <SDL2/SDL_thread.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sndfile.h>
#include <limits.h>

#include <stdint.h>

#define min(A,B) ((A) > (B)? (B) : (A))
#define max(A,B) ((A) > (B)? (A) : (B))

typedef double pixel_t;

static void *malloc_check(size_t len) {
  void *p;
  p = malloc(len);
  if (! p) {
    printf("No mem.\n");
    exit(-1);
  }
  return p;
}

float frandom(void) {
  return (float)(random()) / INT_MAX;
}


#include "images.h"
#include "palettes.h"

#define SEED_VAL (0.5 * PALETTE_LEN)
#define MAX_SEED_R (min_W_H/5)

int n_images = 0;
image_t *images = NULL;

const float minuscule = 1e-3;

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


int W = 1920 / 3;
int H = 1080 / 3;
int min_W_H, max_W_H;
pixel_t *pixbuf = NULL;
int pixbuf_bytes = 0;
pixel_t *apex = NULL;
fftw_complex *pixbuf_f;
fftw_plan plan_backward;
fftw_plan plan_forward;
fftw_complex *apex_f;
fftw_plan plan_apex;

typedef enum {
  ao_left = 2,
  ao_right = 8,
  ao_up = 4,
  ao_down = 1,
} apex_opt_t;

void make_apex(double apex_r, double burn_amount, char apex_opt);

void fft_init(void) {
  fftw_init_threads();
  fftw_plan_with_nthreads(4);

  pixbuf_bytes = W * H * sizeof(pixel_t);

  pixbuf = (pixel_t *) malloc_check(pixbuf_bytes);
  bzero(pixbuf, pixbuf_bytes);

  int half_W = (W / 2) + 1;
  apex = (pixel_t*)malloc_check(pixbuf_bytes);
  bzero(apex, pixbuf_bytes);
  apex_f = fftw_malloc(sizeof(fftw_complex) * H * half_W);
  plan_apex = fftw_plan_dft_r2c_2d(H, W, apex, apex_f, FFTW_ESTIMATE);

  pixbuf_f = fftw_malloc(sizeof(fftw_complex) * H * half_W);
  plan_forward = fftw_plan_dft_r2c_2d(H, W, pixbuf, pixbuf_f, FFTW_ESTIMATE);
  plan_backward = fftw_plan_dft_c2r_2d(H, W, pixbuf_f, pixbuf, FFTW_ESTIMATE);

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

void make_apex(double apex_r, double burn_amount, char apex_opt) {
  int x, y;

  apex_r = min(apex_r, min_W_H/2 - 2);

  double burn_factor = burn_amount + 1.0;

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
  int W2 = W >> 1;
  int H2 = H >> 1;

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
      if (apex_opt) {
        bool neg = (
          ((apex_opt & ao_left) && (x < W2))
          ||
          ((apex_opt & ao_right) && (x > W2))
          ||
          ((apex_opt & ao_up) && (y < H2))
          ||
          ((apex_opt & ao_down) && (y > H2)));
        if (neg) {
          const double vv = -1.8;
          if ((apex_opt & (ao_left | ao_right)) && (apex_opt & (ao_up | ao_down)))
            // two directions pressed simultaneously
            v *= vv;
          else
            v *= vv * 30;
        }
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
    for (x = 0; x < x_fold; x ++) {
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

#define UNPIXELIZE_BITS 5

void render(Uint32 *winbuf, const int winW, const int winH,
            palette_t *palette, pixel_t *pixbuf,
            int multiply_pixels, int colorshift, char pixelize,
            unsigned char invert)
{
  assert((W * multiply_pixels) == winW);
  assert((H * multiply_pixels) == winH);

  int x, y;
  int mx, my;
  int pitch = winW;
  int one_screen_row_pitch = pitch - multiply_pixels;
  int one_multiplied_row_pitch = (pitch - winW)
                                 + (multiply_pixels - 1)*pitch;

  Uint32 *winpos = winbuf;
  pixel_t *pixbufpos = pixbuf;

  int pixelize_mask = ~(INT_MAX << pixelize);
  int pixelize_offset_x = (pixelize_mask - (W & pixelize_mask)) >> 1;
  int pixelize_offset_y = (pixelize_mask - (H & pixelize_mask)) >> 1;

  int invert_mask = INT_MAX;
  invert_mask = ~(invert_mask << UNPIXELIZE_BITS);
  static float invert_offset = 0;
  invert_offset += .014;
  int _invert_offset_x = (int)(((sin(invert_offset) * (winH/2) ) + invert/2));
  int _invert_offset_y = (int)(((cos(invert_offset) * (winW/2) ) + invert/2));

  static float move_pixlz_offset = 0;
  move_pixlz_offset += .1;
  static int pxlz_dir = 0;
  if (move_pixlz_offset > pixelize_mask) {
    move_pixlz_offset = 0;
    pxlz_dir = (pxlz_dir + 1) % 4;
  }

  int _pixelize_offset = move_pixlz_offset * 10;
  pixelize_offset_x = (pixelize_offset_x + _pixelize_offset) & pixelize_mask;
  pixelize_offset_y = (pixelize_offset_y + _pixelize_offset) & pixelize_mask;
  if (pxlz_dir & 1) {
    pixelize_offset_x = -pixelize_offset_x;
  }
  if (pxlz_dir & 2) {
    pixelize_offset_y = -pixelize_offset_y;
  }

#define AVERAGING 0

#if AVERAGING
  pixel_t pmin, pmax, psum;
  pmin = pmax = *pixbufpos;
  psum = 0;
#endif

  for (y = 0; y < H; y++) {
    for (x = 0; x < W; x++, pixbufpos++) {
      pixel_t pix = *pixbufpos;
      if (pix >= palette->len) {
        pix -= palette->len * (int)(pix) / palette->len;
        *pixbufpos = pix;
      }
      else
      if (pix < 0.001) {
        pix = 0;
        //pix += palette->len * (1 + (int)(-pix) / palette->len);
        *pixbufpos = pix;
      }
#if AVERAGING
      if (my == 0) {
        psum += pix;
        pmin = min(pmin, pix);
        pmax = max(pmax, pix);
      }
#endif
      if (pixelize) {
        int xx = (((x + pixelize_offset_x) & ~pixelize_mask) - pixelize_offset_x) + (pixelize_mask >> 1);
        int yy = (((y + pixelize_offset_y) & ~pixelize_mask) - pixelize_offset_y) + (pixelize_mask >> 1);
        pix = *(pixbuf + max(0,min(W-1,xx)) + max(0,min(H-1,yy))*W);
      }

      unsigned int col = (unsigned int)pix + colorshift;
      col %= palette->len;

      Uint32 raw = palette->colors[col];

      if (invert) {
        if ((((x + _invert_offset_x) & invert_mask) <= invert)
            || (((y + _invert_offset_y) & invert_mask) <= invert))
          raw = ~raw;
      }

      Uint32 *p = winpos;

      for (my = 0; my < multiply_pixels; my++) {
        for (mx = 0; mx < multiply_pixels; mx++) {
          *p = raw;
          p ++;
        }
        p += one_screen_row_pitch;
      }

      winpos += multiply_pixels;
    }
    winpos += one_multiplied_row_pitch;
  }

#if AVERAGING
  printf("%.3f %.3f %.3f\r", pmin/PALETTE_LEN, pmax/PALETTE_LEN, psum/((float)W*H*PALETTE_LEN));
  fflush(stdout);
#endif
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

void seed_image(int x, int y, pixel_t *img, int w, int h, pixel_t intensity) {
  pixel_t *pixbuf_pos = pixbuf + y * W + x;
  pixel_t *pixbuf_end = pixbuf + W * H;
  int pixbuf_pitch = max(0, W - w);

  pixel_t *img_pos = img;
  int xx, yy;
  for (yy = 0; yy < h; yy++) {
    for (xx = 0; (xx < w) && (pixbuf_pos < pixbuf_end); xx++) {
      pixel_t add = (*img_pos) * 0.42651 * intensity * PALETTE_LEN;
      (*pixbuf_pos) += add;
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
float want_fps = 25;

Uint32 *winbuf;
SDL_Renderer *renderer;
SDL_Texture *texture;
int winW;
int winH;
palette_t palette;
palette_t blended_palette;
FILE *out_stream = NULL;
FILE *out_params = NULL;
FILE *in_params = NULL;
int in_params_content_start;
int out_params_content_start;
int multiply_pixels = 1;
int colorshift;

typedef struct {
  int random_seed;
  bool start_blank;
} init_params_t;

typedef struct {
  float apex_r;
  char apex_opt;
  float burn_amount;
  float axis_colorshift;
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
  char pixelize;
  char invert;
  int please_drop_img_x;
  int please_drop_img_y;
  float seed_intensity;
  float colorshift_constant;
  float palette_change;
  float min_burn;
} params_t;

const int params_file_id = 0x23315;
const int params_version = 2;

init_params_t ip;
params_t p = {
  .apex_r=3.35,
  .apex_opt = 0,
  .burn_amount = .0013,
  .axis_colorshift = 0.1,
  .do_stop = false,
  .do_go = false,
  .do_wavy = true,
  .do_stutter = false,
  .do_blank = false,
  .do_maximize = false,
  .force_symm = true,
  .symm = symm_x,
  .seed_r = 23,
  .n_seed = 0,
  .wavy_amp = .002,
  .please_drop_img = -1,
  .pixelize = 0,
  .invert = 0,
  .please_drop_img_x = INT_MAX,
  .please_drop_img_y = INT_MAX,
  .seed_intensity = 1,
  .colorshift_constant = 2,
};

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

  float want_frame_period = (want_fps > .1? 1000. / want_fps : 0);
  float last_ticks = (float)SDL_GetTicks() - want_frame_period;

  for (;;) {

    SDL_SemWait(please_render);
    if (! running)
      break;

    while (want_frame_period) {
      int elapsed = SDL_GetTicks() - last_ticks;
      if (elapsed >= want_frame_period) {
        last_ticks += want_frame_period * (int)(elapsed / want_frame_period);
        break;
      }
      SDL_Delay((int)want_frame_period - elapsed);
    }

    if (out_stream) {
      SDL_SemWait(saving_done);
    }

    render(winbuf, winW, winH, &palette, pixbuf, multiply_pixels, colorshift, p.pixelize, p.invert);

    SDL_UpdateTexture(texture, NULL, winbuf, winW * sizeof(Uint32));

    SDL_RenderClear(renderer);
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_RenderPresent(renderer);

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
      fwrite(winbuf, sizeof(Uint32), winW * winH, out_stream);
    }
    SDL_SemPost(saving_done);
  }

  return 0;
}



SNDFILE *audio_sndfile = NULL;
SDL_AudioSpec audio_spec;
int audio_bytes_per_video_frame = 1;
volatile int audio_too = 0;
char *audio_path = NULL;
bool audio_sync_verbose = true;

void audio_play_callback(void *userdata, Uint8 *stream, int len) {
  bool smoothen = false;
  static int audio_bytes_played = 0;
  static short last_frame[2];
  int frame_size = audio_spec.channels * sizeof(short);
  int want_bytes_played = frames_rendered * audio_bytes_per_video_frame;

  // eyes are slower than ears. The sound for the frame should rather come just
  // after the frame is rendered than before.
  want_bytes_played = max(0, want_bytes_played - audio_bytes_per_video_frame);

  int diff = want_bytes_played - audio_bytes_played;

  if (diff < -audio_bytes_per_video_frame) {
    // audio too fast.
    if (audio_sync_verbose)
      printf("audio too fast. %d < %d\n",  diff, -audio_bytes_per_video_frame);
    audio_too = -diff;
    sf_seek(audio_sndfile, diff/frame_size, SEEK_CUR);
    audio_bytes_played = want_bytes_played;
    smoothen = true;
  }
  else
  if (diff > audio_bytes_per_video_frame) {
    // audio too slow.
    if (audio_sync_verbose)
      printf("audio too slow. %d > %d\n", diff, audio_bytes_per_video_frame);
    audio_too = -diff;
    sf_seek(audio_sndfile, diff/frame_size, SEEK_CUR);
    audio_bytes_played = want_bytes_played;
    smoothen = true;
  }

  if (! sf_read_short(audio_sndfile, (short*)stream, len/sizeof(short))) {
    printf("Audio file ended. Stop. (%s)\n", audio_path);
    running = false;
  }
  else
    audio_bytes_played += len;

  if (smoothen) {
    int i;
    int ch;
    short *x = (short*)stream;
    for (i = 0; i < 15; i++) {
      for (ch = 0; ch < audio_spec.channels; ch++) {
        x[(i<<1) + ch] = (((i + 1) * x[(i<<1) + ch]) + ((15 - i) * last_frame[ch])) >> 4;
      }
    }
  }
  memcpy(&last_frame, stream + len - sizeof(last_frame), sizeof(last_frame));
}

int main(int argc, char *argv[])
{
  bool usage = false;
  bool error = false;
  bool cmdline_requests_start_blank = false;
  bool fullscreen = false;

  int c;

  char *out_stream_path = NULL;
  char *out_params_path = NULL;
  char *in_params_path = NULL;

  ip.random_seed = time(NULL);
  ip.start_blank = false;

  int divide_pixels = 1;

  while (1) {
    c = getopt(argc, argv, "bha:d:f:g:m:p:r:u:i:o:O:P:F");
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

      case 'd':
        divide_pixels = atoi(optarg);
        break;

      case 'F':
        fullscreen = true;
        break;

      case 'm':
        multiply_pixels = atoi(optarg);
        break;

      case 'f':
        want_fps = atof(optarg);
        break;

      case 'a':
        p.apex_r = atof(optarg);
        break;

      case 'u':
        p.burn_amount = atof(optarg);
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

      case 'b':
        cmdline_requests_start_blank = true;
        break;

      case 'p':
        audio_path = optarg;
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
"  burnscope -g 320x200 -m 2 -f 30\n"
"\n"
"Options:\n"
"\n"
"  -g WxH   Set animation width and height in number of pixels.\n"
"           Default is '-g %dx%d'.\n"
"  -F       Start in full-screen mode.\n"
"  -f fps   Set desired framerate to <fps> frames per second. The framerate\n"
"           may slew if your system cannot calculate fast enough.\n"
"           If zero, run as fast as possible. Default is %.1f.\n"
"  -m N     Multiply each pixel N times in width and height, to give a larger\n"
"           picture. This will also multiply the window size.\n"
"  -d N     Same as -m, except keeping identical final -g dimensions.\n"
"  -a W     Set apex radius, i.e. the blur distance. Default is %.3f.\n"
"  -u N.n   Set underdampening factor (decimal). Default is %.3f.\n"
"           Reduces normal blur dampening by this factor.\n"
"  -b       Start out blank.\n"
"  -r seed  Supply a random seed to start off with.\n"
"  -O file  Write raw video data to file (grows large quickly). Can be\n"
"           converted to a video file using e.g. ffmpeg.\n"
"  -o file  Write live control parameters to file for later playback, see -i.\n"
"  -i file  Play back previous control parameters (possibly in a different\n"
"           resolution and streaming video to file...)\n"
"  -p file  Play back audio file in sync with actual framerate.\n"
"           The file format should match your sound card output format\n"
"           exactly.\n"
, W, H, want_fps, p.apex_r, p.burn_amount
);
    if (error)
      return 1;
    return 0;
  }

  const int maxpixels = 1e4;

  if ((W < 3) || (W > maxpixels) || (H < 3) || (H > maxpixels)) {
    fprintf(stderr, "width and/or height out of bounds: %dx%d\n", W, H);
    exit(1);
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

  if (((1. + p.burn_amount) > -minuscule) && ((1. + p.burn_amount) < minuscule)) {
    fprintf(stderr, "Underdampening too close to 1.0 (-u). Limit is %f.\n",
        minuscule);
    exit(1);
  }

  winW = W;
  winH = H;

  if (multiply_pixels > 1) {
    winW *= multiply_pixels;
    winH *= multiply_pixels;
  }
  else
    multiply_pixels = 1;

  if (divide_pixels > 1) {
    if ((winW % divide_pixels) || (winH % divide_pixels)) {
      fprintf(stderr, "-d %d: both -g dimensions must be a multiple of -d\n",
              divide_pixels);
      exit(1);
    }
    multiply_pixels = divide_pixels;
    W = winW / multiply_pixels;
    H = winH / multiply_pixels;
  }

  if ( (winW > maxpixels) || (winH > maxpixels) ) {
    fprintf(stderr, "pixel multiplication is too large: %dx%d times %d = %dx%d\n",
            W, H, multiply_pixels, winW, winH);
    exit(1);
  }

  if (out_stream_path) {
    if (access(out_stream_path, F_OK) == 0) {
      fprintf(stderr, "file exists, will not overwrite: %s\n", out_stream_path);
      exit(1);
    }
    out_stream = fopen(out_stream_path, "w");
    audio_sync_verbose = false;
  }

  if (out_params_path) {
    if (access(out_params_path, F_OK) == 0) {
      fprintf(stderr, "file exists, will not overwrite: %s\n", out_params_path);
      exit(1);
    }
    out_params = fopen(out_params_path, "w+");
  }

  if (in_params_path) {
    in_params = fopen(in_params_path, "r");
  }

  printf("burnscope: %dx%d  -->  video: %dx%d\n", W, H, winW, winH);

  int in_params_framelen = sizeof(p);
  int in_params_read_framelen = 0;

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
    int diff;
    diff = ip_len - sizeof(ip);
    if (diff < 0) {
      printf("Parameter file header is smaller: %d instead of %d\n", ip_len, (int)sizeof(ip));
    }
    else
    if (diff > 0) {
      printf("Parameter file header is larger, reading smaller header: %d instead of %d\n", (int)sizeof(ip), ip_len);
      ip_len = sizeof(ip);
    }
    fread(&ip, ip_len, 1, in_params);
    if (diff > 0)
      fseek(in_params, diff, SEEK_CUR);

    fread(&in_params_framelen, sizeof(in_params_framelen), 1, in_params);
    if (in_params_framelen < sizeof(p)) {
      printf("Parameter file frame size is smaller: %d instead of %d\n", in_params_framelen, (int)sizeof(p));
    }
    if (in_params_framelen > sizeof(p)) {
      printf("Parameter file frame size is larger, reading smaller frame: %d instead of %d\n", (int)sizeof(p), in_params_framelen);
    }

    in_params_read_framelen = min(in_params_framelen, sizeof(p));

    in_params_content_start = ftell(in_params);
    printf("in_params start @%d\n", in_params_content_start);
  }

  printf("random seed: %d\n", ip.random_seed);
  srandom(ip.random_seed);

  if (cmdline_requests_start_blank) {
    ip.start_blank = true;
  }

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

    out_params_content_start = ftell(out_params);
  }


  read_images("./images", &images, &n_images, W, H);


  if ( SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_JOYSTICK
                | (audio_path? SDL_INIT_AUDIO : 0))
       < 0 )
  {
    fprintf(stderr, "Unable to init SDL: %s\n", SDL_GetError());
    exit(1);
  }

  SDL_Window *window;
  window = SDL_CreateWindow("burnscope_fft", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                            winW, winH, 0);

  if (!window) {
    fprintf(stderr, "Unable to set %dx%d video: %s\n", winW, winH, SDL_GetError());
    exit(1);
  }

  if (fullscreen)
    SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN);

  renderer = SDL_CreateRenderer(window, -1, 0);
  if (!renderer) {
    fprintf(stderr, "Unable to set %dx%d video: %s\n", winW, winH, SDL_GetError());
    exit(1);
  }

  SDL_ShowCursor(SDL_DISABLE);
  SDL_PixelFormat *pixelformat = SDL_AllocFormat(SDL_PIXELFORMAT_RGBA8888);
  texture = SDL_CreateTexture(renderer, pixelformat->format,
                              SDL_TEXTUREACCESS_STREAMING, winW, winH);
  if (!texture) {
    fprintf(stderr, "Cannot create texture\n");
    exit(1);
  }

  const int n_joysticks = SDL_NumJoysticks();

  printf("%d joysticks were found.\n", n_joysticks);

  SDL_Joystick **joysticks = NULL;

  if (n_joysticks) {
    SDL_JoystickEventState(SDL_ENABLE);

    joysticks = malloc(sizeof(SDL_Joystick*) * n_joysticks);

    int i;
    for (i = 0; i < n_joysticks; i++)
    {
      printf("%2d: '%s'\n", i, SDL_JoystickNameForIndex(i));

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

  make_palettes(pixelformat);
  make_palette(&palette, PALETTE_LEN,
               palette_defs[0],
               pixelformat);

  make_palette(&blended_palette, PALETTE_LEN,
               palette_defs[0],
               pixelformat);

  fft_init();

  winbuf = malloc_check(winW * winH * sizeof(Uint32));

  if (! ip.start_blank) {
    int i, j;
    j = 2*p.apex_r + 1;
    j *= j;
    j = W * H / j;
    for (i = 0; i < j; i ++) {
      seed(pixbuf, W, H, random() % (W), random() % (H), SEED_VAL, p.apex_r);
    }
  }
  else {
    printf("blank\n");
  }

  float wavy = 0;
  bool do_print = true;
  float wavy_speed = .5;
  double use_burn = 1.002;

  please_render = SDL_CreateSemaphore(0);
  please_save = SDL_CreateSemaphore(0);
  rendering_done = SDL_CreateSemaphore(0);
  saving_done = SDL_CreateSemaphore(1);

  SDL_Thread *render_thread_token = SDL_CreateThread(render_thread, NULL, "render");
  SDL_Thread *save_thread_token = NULL;
  if (out_stream)
    SDL_CreateThread(save_thread, NULL, "save");

  fftw_execute(plan_forward);

  if (audio_path) {
    SF_INFO audio_sndfile_info;
    audio_sndfile = sf_open(audio_path, SFM_READ, &audio_sndfile_info);

    if (audio_sndfile == NULL) {
      fprintf(stderr, "Cannot open audio file: %s\n", audio_path);
      exit(1);
    }

    int frames = audio_sndfile_info.frames;
    int seconds = frames / audio_sndfile_info.samplerate;
    int minutes = seconds / 60;
    seconds %= 60;
    printf("%s:\n  %d ch %d Hz %dm%ds",
           audio_path,
           audio_sndfile_info.channels,
           audio_sndfile_info.samplerate,
           minutes, seconds
        );
    if (audio_sndfile_info.format == SF_FORMAT_PCM_16) {
      printf(" 16bit WAV");
    }
    printf("\n");

    SDL_AudioSpec audio_want_spec;
    audio_want_spec.freq = audio_sndfile_info.samplerate;
    audio_want_spec.format = AUDIO_S16;
    audio_want_spec.channels = audio_sndfile_info.channels;
    audio_want_spec.samples = audio_want_spec.freq / (want_fps > .1? want_fps : 25);
    audio_want_spec.callback = audio_play_callback;
    audio_want_spec.userdata = NULL;

    if (SDL_OpenAudio(&audio_want_spec, &audio_spec) < 0) {
      fprintf(stderr, "Cannot open audio: %s\n", SDL_GetError());
      exit(1);
    }

    printf("Opened audio:\n"
           "  %d ch %d Hz %d buf",
           audio_spec.channels,
           audio_spec.freq,
           audio_spec.samples);
    if (audio_spec.format == AUDIO_S16) {
      printf(" 16bit WAV");
    }
    printf("\n");
    audio_bytes_per_video_frame = audio_spec.freq * audio_spec.channels * sizeof(short)
                                  / (want_fps > .1? want_fps : 25);
    printf("audio: %d bytes per video frame.\n", audio_bytes_per_video_frame);
    SDL_PauseAudio(0);
  }

  double apex_r_center = p.apex_r;
  double burn_center = p.burn_amount;
  double burn_jump = 0;
  bool do_back = false;
  bool do_rerecord_params = false;
  bool do_rerecord_params_overlay = false;
  bool do_rerecord_params_overlay_started = false;
  int had_outparams = 0;
  int img_seeding = -1;
  int img_seeding_slew = 0;

  while (running)
  {
#define BACK_SPEED 4
#define BACK_SEEK (BACK_SPEED + 1)
    if (do_back && (frames_rendered > (BACK_SEEK+1))) {
      if (in_params) {
        fseek(in_params, -BACK_SEEK * in_params_framelen, SEEK_CUR);
      }
      if (out_params) {
        fseek(out_params, -BACK_SEEK * sizeof(p), SEEK_CUR);
      }
      frames_rendered -= BACK_SEEK;
    }

    if (in_params || out_params) {
      // read parameters from in_params (FILE*).
      // If out_params (FILE*) is being written and the user has rewound back
      // to a place where parameters are already written to out_params, read
      // those.
      // If overlay is requested, the live joystick commands are partly kept
      // instead of being overwritten by the parameters read from file. In that
      // case read to a separate location first, and only overwrite those
      // values that have wiggled on the joystick since overlay recording
      // started.
      static params_t p_from_file;
      static params_t prev_p;
      static params_t overlay_mask;

      int in_params_skip_bytes = 0;
      params_t *read_in_params_to = NULL;

      if (do_rerecord_params && do_rerecord_params_overlay) {
        if (do_rerecord_params_overlay_started) {
          bzero(&overlay_mask, sizeof(overlay_mask));
          do_rerecord_params_overlay_started = false;
        }
        else {
          #define check_param(name) \
            if (prev_p.name != p.name) {overlay_mask.name = 1; printf("%s\n", #name); }

          check_param(do_stop);
          check_param(apex_r);
          check_param(apex_opt);
          check_param(burn_amount);
          check_param(axis_colorshift);
          check_param(do_stop);
          check_param(do_go);
          check_param(do_wavy);
          check_param(do_stutter);
          check_param(do_blank);
          check_param(do_maximize);
          check_param(force_symm);
          check_param(symm);
          check_param(seed_r);
          check_param(n_seed);
          check_param(wavy_amp);
          check_param(please_drop_img);
          check_param(pixelize);
          check_param(invert);
          check_param(please_drop_img_x);
          check_param(please_drop_img_y);
          check_param(seed_intensity);
          check_param(colorshift_constant);
          check_param(palette_change);

          #undef check_param
        }
      }

      if (do_rerecord_params) {
        if (do_rerecord_params_overlay)
          read_in_params_to = &p_from_file;
        else
          in_params_skip_bytes = in_params_framelen;
      }
      else
        read_in_params_to = &p;

      if (read_in_params_to) {
        bool do_read_in_params = true;

        if (out_params && (had_outparams > frames_rendered)) {
          // user has previously played these params and recorded them to the
          // out_params file, and has possibly recorded new parameters in the
          // process. Read those that were written earlier.
          if (fread(read_in_params_to, sizeof(p), 1, out_params)) {
            fseek(out_params, -sizeof(p), SEEK_CUR);
            do_read_in_params = false;
          }
        }
        if (! do_read_in_params) {
          in_params_skip_bytes = in_params_framelen;
        }
        else
        if (in_params) {
          if (! fread(read_in_params_to, in_params_read_framelen, 1, in_params)) {
            printf("End of input parameters. Stop. (%s)\n", in_params_path);
            running = false;
          }
          else
          if (in_params_framelen > in_params_read_framelen) {
            // skip trailing bytes if params file frame is larger than my params_t.
            in_params_skip_bytes = in_params_framelen - in_params_read_framelen;
          }
        }
      }

      if (in_params_skip_bytes && in_params)
        fseek(in_params, in_params_skip_bytes, SEEK_CUR);

      if (do_rerecord_params && do_rerecord_params_overlay) {
        #define check_param(name) \
          if (! overlay_mask.name) p.name = p_from_file.name

        check_param(do_stop);
        check_param(apex_r);
        check_param(apex_opt);
        check_param(burn_amount);
        check_param(axis_colorshift);
        check_param(do_stop);
        check_param(do_go);
        check_param(do_wavy);
        check_param(do_stutter);
        check_param(do_blank);
        check_param(do_maximize);
        check_param(force_symm);
        check_param(symm);
        check_param(seed_r);
        check_param(n_seed);
        check_param(wavy_amp);
        check_param(please_drop_img);
        check_param(pixelize);
        check_param(invert);
        check_param(please_drop_img_x);
        check_param(please_drop_img_y);
        check_param(seed_intensity);
        check_param(colorshift_constant);
        check_param(palette_change);

        #undef check_param

        memcpy(&prev_p, &p, sizeof(p));
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

    colorshift = normalize_colorshift;

    {
      static float colorshift_accum = 0;
      colorshift_accum += ((float)PALETTE_LEN / 120) * p.axis_colorshift * fabs(p.axis_colorshift);
      colorshift_accum += p.colorshift_constant;
      colorshift += (int)colorshift_accum;
    }
    colorshift %= PALETTE_LEN;

    {
      float palette_change_per_frame = (0.002 + p.palette_change) / (want_fps > .1? want_fps : 25);
      static int left_pal_idx = 0;
      static int right_pal_idx = 1;
      static float palette_blend_position = 0;

      palette_blend_position += palette_change_per_frame;

      if (palette_blend_position > 1.) {
        left_pal_idx = right_pal_idx;
        right_pal_idx = (left_pal_idx + 1) % n_palettes;
        palette_blend_position -= 1;
      }
      else
      if (palette_blend_position < 0) {
        right_pal_idx = left_pal_idx;
        left_pal_idx = (right_pal_idx > 0? right_pal_idx - 1 : n_palettes - 1);
        palette_blend_position += 1;
      }

      palette_t *left_pal = &palettes[left_pal_idx % n_palettes];
      palette_t *right_pal = &palettes[right_pal_idx % n_palettes];

      blend_palettes(&palette, left_pal, right_pal, palette_blend_position);
    }

    if (img_seeding >= 0) {
      if (img_seeding_slew) {
        img_seeding_slew --;
      }
      else {
        int _img_seeding = img_seeding;
        if (_img_seeding >= 10) {
          _img_seeding -= 10;
          if (_img_seeding < n_images) {
            p.please_drop_img_x = -images[_img_seeding].width / 2;
            p.please_drop_img_y = -images[_img_seeding].height / 2;
          }
        }
        p.please_drop_img = _img_seeding;
        img_seeding_slew = 0; // 2
      }
    }


    if (do_calc) {
      float t = (float)frames_rendered / 100.;
      wavy = sin(wavy_speed*t);

      use_burn = p.burn_amount;

      if (p.do_wavy) {
        use_burn += (p.wavy_amp * wavy);
      }

      if ((p.burn_amount > 0) && (p.apex_r > 10))
        use_burn += (p.apex_r - 10) * .0000625;
    }

    if (out_params) {
      fwrite(&p, sizeof(p), 1, out_params);
      had_outparams = max(had_outparams, frames_rendered);
    }

    // clear those values that were handled above but still needed to be
    // recorded to out_params:
    {
      if (p.do_blank) {
        p.do_blank = false;
      }
      if (p.do_maximize) {
        p.do_maximize = false;
      }
    }

    if (do_calc) {
      p.n_seed = max(0, min(100, p.n_seed));
      while (running && p.n_seed) {
        p.n_seed --;
        int seedx = random() % W;
        int seedy = random() % H;
        seed(pixbuf, W, H, seedx, seedy, SEED_VAL, p.seed_r);

        if ((p.symm == symm_x) || (p.symm == symm_xy))
          // seedx = 0 ==> seedx = W -1
          seed(pixbuf, W, H, W-1 - seedx, seedy, SEED_VAL, p.seed_r);

        if ((p.symm == symm_y) || (p.symm == symm_xy))
          seed(pixbuf, W, H, seedx, H-1 - seedy, SEED_VAL, p.seed_r);
        if (p.symm == symm_point)
          seed(pixbuf, W, H, W-1 - seedx, H-1 - seedy, SEED_VAL, p.seed_r);
      }

      if (p.please_drop_img >= 0) {
        if (p.please_drop_img < n_images) {
          int W2 = W >> 1;
          int H2 = H >> 1;
          image_t *img = &images[p.please_drop_img];

          if (p.please_drop_img_x == INT_MAX)
            p.please_drop_img_x = random() % (30 + W - img->width);
          else
            p.please_drop_img_x += W2;

          if (p.please_drop_img_y == INT_MAX)
            p.please_drop_img_y = random() % (30 + H - img->height);
          else
            p.please_drop_img_y += H2;

          float intensity = p.seed_intensity;
          intensity = .001 + .1 * intensity * intensity;
#if 0
          printf("img x%d y%d %d %f\n", p.please_drop_img_x, p.please_drop_img_y,
                 p.please_drop_img, intensity);
#endif

          seed_image(p.please_drop_img_x, p.please_drop_img_y, img->data, img->width, img->height,
                     1.);
        }
        p.please_drop_img = -1;
        p.please_drop_img_x = INT_MAX;
        p.please_drop_img_y = INT_MAX;
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
        pixel_t *pf = pixbuf_f[x];
        pixel_t *af = apex_f[x];
        pixel_t a, b, c, d;
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

        if (p.symm != symm_none)
          p.force_symm = true;
      }
    }

    static char printcount = 0;
    if (printcount++ >= 10) {
      if (printcount >= 50)
        do_print = true;

      if (do_print) {
        do_print = false;
        printcount = 0;
        printf("%.1ffps apex_r=%f_opt%d burn=%f(%f) audio_sync=%d\n",
               1000./(avg_frame_period>>AVG_SHIFTING),
               p.apex_r,p.apex_opt,
               use_burn,
               p.burn_amount,
               audio_too);
        audio_too = 0;
        fflush(stdout);
      }
    }

    while (running) {
      SDL_Event event;
      while (SDL_PollEvent(&event))
      {

        switch (event.type)
        {
          case SDL_KEYDOWN:
            {
              int c = event.key.keysym.sym;

              switch(c) {
                case SDLK_ESCAPE:
                  printf("Escape key. Stop.\n");
                  running = false;
                  break;

                case SDLK_RIGHT:
                  p.burn_amount += .0002;
                  break;
                case SDLK_LEFT:
                  p.burn_amount -= .0002;
                  break;
                case SDLK_UP:
                  p.wavy_amp += .0001;
                  break;
                case SDLK_DOWN:
                  p.wavy_amp -= .0001;
                  break;

                case ' ':
                  p.n_seed ++;
                  break;

                case 8: // del
                  p.do_blank = true;
                  break;

                case 'm':
                  p.symm = (p.symm + 1) % SYMMETRY_KINDS;
                  p.force_symm = true;
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
                  p.burn_amount -= .002;
                  break;
                case 'w':
                  p.burn_amount -= .0003;
                  break;
                case 'e':
                  p.burn_amount = .005;
                  p.apex_r = 8.01 * min_W_H / 240.;
                  break;
                case 'r':
                  p.burn_amount += .0003;
                  break;
                case 't':
                  p.burn_amount += .002;
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

                case 'l':
                  wavy_speed += .5;
                  break;

                case 'k':
                  wavy_speed -= .5;
                  break;

                case 'y':
                  img_seeding = 14;
                  break;

                case 'u':
                  img_seeding = 10;
                  break;

                case 'i':
                  img_seeding = 11;
                  break;

                case 'o':
                  img_seeding = 12;
                  break;

                case 'p':
                  img_seeding = 13;
                  break;

                case ']':
                  p.axis_colorshift += 0.1;
                  break;

                case '[':
                  p.axis_colorshift -= 0.1;
                  break;

                default:

                  if ((c >= '1') && (c <= '9')) {
                    int img_idx = c - '1';
                    img_seeding = img_idx;
                  }
                  else
                    printf("keysym = %c (%d)\n", (char)max(0x20,c), c);
                  break;
              }


            }
            do_print = true;
            break;

          case SDL_KEYUP:
            {
              int c = event.key.keysym.sym;

              int img_seeding_off = -2;
              switch (c) {
                case 'u':
                  img_seeding_off = 10;
                  break;

                case 'i':
                  img_seeding_off = 11;
                  break;

                case 'o':
                  img_seeding_off = 12;
                  break;

                case 'p':
                  img_seeding_off = 13;
                  break;

                case 'y':
                  img_seeding_off = 14;
                  break;

                case 13:
                  fullscreen = !fullscreen;
                  if (fullscreen)
                    SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN);
                  else
                    SDL_SetWindowFullscreen(window, 0);
                  break;

                default:
                  if ((c >= '1') && (c <= '9')) {
                    img_seeding_off = c - '1';
                  }
                  break;
              }

              if (img_seeding == img_seeding_off) {
                img_seeding = -1;
                img_seeding_slew = 0;
              }
            }
            break;

          case SDL_QUIT:
            printf("SDL_QUIT. Stop.\n");
            running = false;
            break;


#define DO_JOY_DEBUG 0

#if DO_JOY_DEBUG
#define JOYMSG(fmt, args...) printf(fmt, ##args )
#else
#define JOYMSG(fmt, args...)
#endif

          case SDL_JOYAXISMOTION:
            {
              float axis_val = event.jaxis.value;
              axis_val /= 32768;
              // JOYMSG("%d axis %d val %f\n", event.jaxis.which, event.jaxis.axis, axis_val);

#define calc_axis_val(start, center, end, axis_val) \
                  ((axis_val <= 0)? \
                    ((start) + ((center) - (start)) * (1 + (axis_val))) \
                    : \
                    ((center) + ((end) - (center)) * (axis_val) * (axis_val)))
#define AXIS_MIN .05

                  switch(event.jaxis.axis)
                  {
                    case 3:
                      p.axis_colorshift = -axis_val;
                      break;

                    case 0:
                      p.palette_change = (fabs(axis_val) > .1)? axis_val : 0;
                      break;


                    case 4:
                      {
                        double apex_r_min = max(2.03, apex_r_center * 0.1);
                        double apex_r_max = min(min_W_H/6, apex_r_center * 5);
                        p.apex_r = calc_axis_val(apex_r_min, apex_r_center, apex_r_max, axis_val);
                        if (axis_val < -.5) {
                          JOYMSG("val %f -> maximize\n", axis_val);
                          p.do_maximize = true;
                        }
                      }
                      break;

                    case 1:
                      {
                        const double bmin = -.007;
                        const double bmax = .025;
                        double _burn_center = burn_center + burn_jump;
                        double burn_min = bmin + (_burn_center - bmin) * .2;
                        double burn_max = min(bmax, bmin + ((_burn_center - bmin) * 2));
                        axis_val *= axis_val * axis_val;
                        p.burn_amount = calc_axis_val(burn_min, _burn_center, burn_max, (-axis_val));
                        do_print = true;
                        if (axis_val < -.5) {
                          JOYMSG("val %f -> maximize\n", axis_val);
                          p.do_maximize = true;
                        }
                      }
                      break;


                    case 5:
                      p.invert = ((axis_val + 1) / 2) * ((1 << UNPIXELIZE_BITS) ); \
                      break;

                    case 2:
                      p.pixelize = ((axis_val + 1) / 2) * (min_W_H > 500? 9 : 6);
                      break;

                    default:
                      printf("axis %d = %.3f\n", event.jaxis.axis, axis_val);
                      break;
                  }

              }
            break; // SDL_JOYAXISMOTION

          case SDL_JOYBUTTONDOWN:
              JOYMSG("%d button %d down\n", event.jbutton.which, event.jbutton.button);
              {

                switch (event.jbutton.button) {
                case 3:
                case 0:
                  p.symm = (p.symm + 1) % SYMMETRY_KINDS;
                  p.force_symm = true;
                  p.n_seed = 1;
                  JOYMSG("symmetry: %s\n", symmetry_name[p.symm % SYMMETRY_KINDS]);
                  break;

                case 1:
                  {
                    static int next_image = -1;
                    next_image = (next_image + 1) % n_images;
                    p.please_drop_img = next_image;
                    p.symm = symm_none;
                  }
                  break;

                case 2:
                  p.n_seed = 1;
                  break;

                case 4:
                  p.do_maximize = true;
                  break;

                case 5:
                  p.do_blank = true;
                  break;

                case 6:
                  do_back = true;
                  do_rerecord_params = false;
                  break;

                case 7:
                  do_rerecord_params = ! do_rerecord_params;
                  do_rerecord_params_overlay = do_rerecord_params;
                  do_rerecord_params_overlay_started = do_rerecord_params_overlay;
                  break;

                case 8:
                  do_rerecord_params = ! do_rerecord_params;
                  do_rerecord_params_overlay = false;
                  break;

                default:
                  printf("%2d: button %d = %s\n",
                         event.jbutton.which, event.jbutton.button,
                         event.jbutton.state == SDL_PRESSED ? "pressed" : "released");
                  break;
                }
            }
            break;


          case SDL_JOYBUTTONUP:
            {
              switch (event.jbutton.button) {

              default:
                printf("%2d: button %d = %s\n",
                       event.jbutton.which, event.jbutton.button,
                       event.jbutton.state == SDL_PRESSED ? "pressed" : "released");
                break;
              }
            }
            break;

          case SDL_JOYHATMOTION:  /* Handle Hat Motion */
            p.apex_opt = event.jhat.value;
            break;

          case SDL_JOYBALLMOTION:  /* Handle Joyball Motion */
            printf("%2d: hat %d += %d, %d\n",
                   event.jball.which, event.jball.ball,
                   event.jball.xrel, event.jball.yrel);
            break;

        } // switch event.type, keyboard and sys events

        {
          static bool was_do_rerecord_params = false;
          if (was_do_rerecord_params != do_rerecord_params) {
            if (do_rerecord_params) {
              if (do_rerecord_params_overlay)
                printf("\n o  RECORD OVERLAY  o  o  o  o  o  o  o  o  o  o  o  o\n\n");
              else
                printf("\n o  RECORD    o  o  o  o  o  o  o  o  o  o  o  o  o  o\n\n");
            }
            else
              printf("\n >  playback  >  >  >  >  >  >  >  >  >  >  >  >  >  >\n\n");
            was_do_rerecord_params = do_rerecord_params;
          }
        }

        if (! running)
          break;

      } // while sdl poll event

      if (SDL_SemTryWait(rendering_done) == 0)
        break;
      else
        SDL_Delay(5);
    } // while running, for event polling / idle waiting

  } // while running

  if (running)
    printf("Main loop exited. Stop.\n");
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
        "ffmpeg -vcodec rawvideo -f rawvideo -pix_fmt rgb32 -s %dx%d -i %s ",
        winW, winH, out_stream_path);
    if (audio_path)
      printf("-i %s -acodec ac3 ", audio_path);
    printf("-vcodec mpeg4 -q 1 %s.%d.mp4\n", out_stream_path, winH);
  }
  if (out_params) {
    if (in_params) {
      printf("Copying remaining parameters stream from in to out: %s --> %s\n",
             in_params_path, out_params_path);

      fseek(out_params, 0, SEEK_END);
      int out_has_frames = (ftell(out_params) - out_params_content_start)
                           / sizeof(p);
      fseek(in_params, in_params_content_start
                       + out_has_frames * in_params_framelen, SEEK_SET);

      int copied = 0;
      while (true) {
        if (! fread(&p, in_params_read_framelen, 1, in_params)) {
          break;
        }
        if (! fwrite(&p, sizeof(p), 1, out_params)) {
          printf("Can't write to %s\n", out_params_path);
          break;
        }
        copied ++;
        if (in_params_framelen > in_params_read_framelen) {
          // skip trailing bytes if params file frame is larger than my params_t.
          fseek(in_params, in_params_framelen - in_params_read_framelen, SEEK_CUR);
        }
      }
      printf("copied %d frames.\n", copied);
    }

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

// vim: ts=2 sw=2 et
