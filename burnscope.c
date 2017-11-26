/* burnscope.c
 * (c) 2014 Neels Hofmeyr <neels@hofmeyr.de>
 *
 * This file is part of burnscope, published under the GNU General Public
 * License v3.
 */

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <png.h>
#include <stdint.h>

#define PALETTE_LEN_BITS 12
#define PALETTE_LEN (1 << PALETTE_LEN_BITS)

typedef Uint32 pixel_t;

typedef struct {
  Uint32 *colors;
  unsigned int len;
  SDL_PixelFormat *format;
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

void set_color(palette_t *palette, int i, float r, float g, float b) {
  if (i >= palette->len)
    return;
  palette->colors[i] = SDL_MapRGB(palette->format, r * 255, g * 255, b * 255);
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
  palette->format = format;


  if (n_points < 1) {
    for (i = 0; i < palette->len; i++) {
      float val = (float)i / palette->len;
      set_color(palette, i, val, val, val);
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

      set_color(palette, color_pos, r, g, b);
    }

    p = *next_p;
  }
}


#define min(A,B) ((A) > (B)? (B) : (A))
#define max(A,B) ((A) > (B)? (A) : (B))

#define SUM_RANGE_BITS 8

Uint32 rectangle_sum(pixel_t *pixbuf, int W, int H,
                     int x_start, int y_start,
                     int x_end, int y_end,
                     bool wrap_borders) {
  Uint32 sum = 0;

  if (x_start < 0) {
    if (wrap_borders)
      sum += rectangle_sum(pixbuf, W, H, W - (-x_start), y_start, W, y_end, wrap_borders);
    x_start = 0;
  }
  if (x_end > W) {
    if (wrap_borders)
      sum += rectangle_sum(pixbuf, W, H, 0, y_start, x_end - W, y_end, wrap_borders);
    x_end = W;
  }

  if (y_start < 0) {
    if (wrap_borders)
      sum += rectangle_sum(pixbuf, W, H, x_start, H - (-y_start), x_end, H, wrap_borders);
    y_start = 0;
  }
  if (y_end > H) {
    if (wrap_borders)
      sum += rectangle_sum(pixbuf, W, H, x_start, 0, x_end, y_end - H, wrap_borders);
    y_end = H;
  }

  pixel_t *bufpos = &pixbuf[x_start + W*y_start];
  int pitch = W - (x_end - x_start);
  int xpos, ypos;
  for (ypos = y_start; ypos < y_end; ypos ++) {
    for (xpos = x_start; xpos < x_end; xpos ++) {
      sum += (*bufpos) >> SUM_RANGE_BITS;
      bufpos ++;
    }
    bufpos += pitch;
  }
  return sum;
}

Uint32 surrounding_sum(pixel_t *pixbuf, const int W, const int H,
                       const int x, const int y, const int apex_r,
                       bool wrap_borders) {
  int x_start = x - apex_r;
  int y_start = y - apex_r;
  int wh = 2 * apex_r + 1;
  int x_end = x_start + wh;
  int y_end = y_start + wh;
  return rectangle_sum(pixbuf, W, H, x_start, y_start, x_end, y_end,
                       wrap_borders);
}

void burn(pixel_t *srcbuf, pixel_t *destbuf, const int W, const int H,
          const int apex_r, float divider, const int palette_len,
          bool wrap_borders,
          int rect_x, int rect_y, int rect_w, int rect_h) {
  int x, y;
  int x_end = rect_x + rect_w;
  int y_end = rect_y + rect_h;
  Uint64 sum;
  int pitch = W - rect_w;
  pixel_t *destpos = destbuf + (rect_y * W) + rect_x;
  for (y = rect_y; y < y_end; y++) {
    for (x = rect_x; x < x_end; x++) {
      sum = surrounding_sum(srcbuf, W, H, x, y, apex_r, wrap_borders) / divider;

      *(destpos ++) = sum << SUM_RANGE_BITS;
    }
    destpos += pitch;
  }
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


void render(Uint32 *winbuf, const int winW, const int winH,
            palette_t *palette, pixel_t *pixbuf, const int W, const int H,
            int multiply_pixels, int colorshift, FILE *out_stream, const char *png_out_path)
{
  assert((W * multiply_pixels) == winW);
  assert((H * multiply_pixels) == winH);

  int x, y;
  int mx, my;

  FILE * png_fp;
  png_structp png_ptr = NULL;
  png_infop png_info = NULL;
  png_byte ** png_row_pointers = NULL;
  /* The following number is set by trial and error only. I cannot
     see where it it is documented in the libpng manual.
  */
  const int png_pixel_size = 3;
  const int png_depth = 8;

  if (png_out_path) {
    png_fp = fopen (png_out_path, "wb");
    if (! png_fp) {
      fprintf(stderr, "Cannot open for writing: '%s'", png_out_path);
      return;
    }

    png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
        goto png_create_write_struct_failed;
    }

    png_info = png_create_info_struct (png_ptr);
    if (png_info == NULL) {
        goto png_create_info_struct_failed;
    }

    /* Set up error handling. */

    if (setjmp (png_jmpbuf (png_ptr))) {
        goto png_failure;
    }

    /* Set image attributes. */

    png_set_IHDR (png_ptr,
                  png_info,
                  winW,
                  winH,
                  png_depth,
                  PNG_COLOR_TYPE_RGB,
                  PNG_INTERLACE_NONE,
                  PNG_COMPRESSION_TYPE_DEFAULT,
                  PNG_FILTER_TYPE_DEFAULT);

    png_row_pointers = png_malloc (png_ptr, winH * sizeof (png_byte *));
  }

  Uint32 *winpos = winbuf;
  pixel_t *pixbufpos = pixbuf;
  int winy = 0;
  for (y = 0; y < H; y++) {
    pixel_t *pixbuf_y_pos = pixbufpos;

    for (my = 0; my < multiply_pixels; my ++, winy++) {
      png_byte *png_row;
      if (png_out_path) {
        png_row = png_malloc (png_ptr, sizeof (uint8_t) * winW * png_pixel_size);
        png_row_pointers[winy] = png_row;
      }

      pixbufpos = pixbuf_y_pos;
      int winx = 0;
      for (x = 0; x < W; x++) {
        Uint32 col = (*pixbufpos) + (colorshift << (32 - PALETTE_LEN_BITS));
        col >>= 32 - PALETTE_LEN_BITS;
        Uint32 raw = palette->colors[col];
        Uint8 r, g, b;

        if (png_out_path)
          SDL_GetRGB(raw, palette->format, &r, &g, &b);

        for (mx = 0; mx < multiply_pixels; mx++, winx++) {
          *winpos = raw;
          winpos ++;

          if (png_out_path) {
            *png_row++ = r;
            *png_row++ = g;
            *png_row++ = b;
          }
        }
        pixbufpos ++;
      }
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
      ((Uint32*)winbuf)[i] = palette->colors[i];
    }
  }
#endif

  if (out_stream) {
    fwrite(winbuf, sizeof(pixel_t), winW * winH, out_stream);
  }

  if (png_out_path) {
    png_init_io (png_ptr, png_fp);
    png_set_rows (png_ptr, png_info, png_row_pointers);
    png_write_png (png_ptr, png_info, PNG_TRANSFORM_IDENTITY, NULL);

    for (y = 0; y < winH; y++) {
        png_free (png_ptr, png_row_pointers[y]);
    }
    png_free (png_ptr, png_row_pointers);

 png_failure:
 png_create_info_struct_failed:
    png_destroy_write_struct (&png_ptr, &png_info);
 png_create_write_struct_failed:
    fclose (png_fp);
  }
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


int main(int argc, char *argv[])
{
  int W = 320;
  int H = 240;
  int multiply_pixels = 1;
  int apex_r = 2;
  float underdampen = .996;
  int frame_period = 70;
  bool usage = false;
  bool error = false;
  bool wrap_borders = true;
  bool start_blank = false;
  symmetry_t symm = symm_x;

  int c;
  int random_seed = time(NULL);

  char *out_stream_path = NULL;
  FILE *out_stream = NULL;
  char *png_out_dir = NULL;

  while (1) {
    c = getopt(argc, argv, "a:g:m:p:r:u:O:P:AbBh");
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
        apex_r = atoi(optarg);
        break;

      case 'u':
        underdampen = atof(optarg);
        break;

      case 'r':
        random_seed = atoi(optarg);
        break;

      case 'O':
        out_stream_path = optarg;
        break;

      case 'P':
        png_out_dir = optarg;
        break;

      case 'b':
        wrap_borders = false;
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
"  -a W     Set apex radius, i.e. the blur distance. Default is %d.\n"
"  -u N.n   Set underdampening factor (decimal). Default is %.3f.\n"
"           Reduces normal blur dampening by this factor.\n"
"  -r seed  Supply a random seed to start off with.\n"
"  -b       Assume zeros around borders. Default is to wrap around borders.\n"
"  -B       Start out blank. (Use 's' key to plant seeds while running.)\n"
, frame_period, apex_r, underdampen
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
    int was_apex_r = apex_r;
    int max_dim = max(W, H);
    apex_r = min(max_dim, apex_r);
    apex_r = max(1, apex_r);
    if (apex_r != was_apex_r) {
      fprintf(stderr, "Invalid apex radius (-a). Forcing %d.", apex_r);
    }
  }

  float minuscule = 1e-3;
  if ((underdampen > -minuscule) && (underdampen < minuscule)) {
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

  if (png_out_dir) {
    int err = mkdir(png_out_dir, 0777);
    if (err) {
      fprintf(stderr, "cannot make dir: %s\n", png_out_dir);
    }
  }


  if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) < 0) {
    fprintf(stderr, "Unable to init SDL: %s\n", SDL_GetError());
    exit(1);
  }

  SDL_Window *window;
  window = SDL_CreateWindow("burnscope", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                            winW, winH, 0);

  if (!window) {
    fprintf(stderr, "Unable to set %dx%d video: %s\n", winW, winH, SDL_GetError());
    exit(1);
  }

  SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, 0);
  if (!renderer) {
    fprintf(stderr, "Unable to set %dx%d video: %s\n", winW, winH, SDL_GetError());
    exit(1);
  }

  SDL_ShowCursor(SDL_DISABLE);
  SDL_PixelFormat *pixelformat = SDL_AllocFormat(SDL_PIXELFORMAT_RGBA8888);
  SDL_Texture *texture = SDL_CreateTexture(renderer, pixelformat->format,
                                           SDL_TEXTUREACCESS_STREAMING, winW, winH);
  if (!texture) {
    fprintf(stderr, "Cannot create texture\n");
    exit(1);
  }

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
               pixelformat);

  pixel_t *buf1 = malloc_check(W * H * sizeof(pixel_t));
  pixel_t *buf2 = malloc_check(W * H * sizeof(pixel_t));
  Uint32 *winbuf = malloc_check(winW * winH * sizeof(Uint32));
  bzero(buf1, W * H * sizeof(pixel_t));
  bzero(buf2, W * H * sizeof(pixel_t));

  pixel_t *pixbuf = buf1;
  pixel_t *swapbuf = buf2;

  printf("random seed: %d\n", random_seed);
  srandom(random_seed);


  if (! start_blank) {
    int i, j;
    j = 2*apex_r + 1;
    j *= j;
    j = W * H / j;
    for (i = 0; i < j; i ++) {
      seed(pixbuf, W, H, random() % (W), random() % (H), 0x80000000, apex_r);
    }
  }

  int last_ticks = SDL_GetTicks() - frame_period;
  bool seed_key_down = false;
  int do_seed = 0;

  int frames_rendered = 0;
  int cc = 0;
  float wavy = 0;
  float wavy_amp = .006;
  int colorshift = 0;
  bool running = true;

  while (running)
  {
    bool do_render = false;

    float t = (float)frames_rendered / 100.;
    wavy = sin(t);
    colorshift = palette.len * (0.5 + 0.5 * cos(t*M_PI/50));

    float dampen = underdampen + (wavy_amp * wavy);
    if ((++cc) > 40) {
      //printf("%.5f + %.5f*%.1f = %.5f  apex_r=%3d  colorshift=%6d/%6d    \r", underdampen, wavy_amp, wavy, dampen, apex_r, colorshift, palette.len);
      //fflush(stdout);
      cc = 0;
    }

    float divider = 1 + 2 * apex_r;
    divider *= divider * dampen;

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
        seed(pixbuf, W, H, seedx, seedy, 0x80000000, apex_r);
        if ((symm == symm_x) || (symm == symm_xy))
          seed(pixbuf, W, H, W - seedx, seedy, 0x80000000, apex_r);
        if ((symm == symm_y) || (symm == symm_xy))
          seed(pixbuf, W, H, seedx, H - seedy, 0x80000000, apex_r);
        if (symm == symm_point)
          seed(pixbuf, W, H, W - seedx, H - seedy, 0x80000000, apex_r);
      }

      pixel_t *tmp = swapbuf;
      swapbuf = pixbuf;
      pixbuf = tmp;

      int ww = W;
      int hh = H;
      if ((symm == symm_x) || (symm == symm_xy))
        ww = W - (W >> 1);
      if ((symm == symm_y) || (symm == symm_xy) || (symm == symm_point))
        hh = H - (H >> 1);

      burn(swapbuf, pixbuf, W, H, apex_r, divider, palette.len, wrap_borders,
           0, 0, ww, hh);

      if (symm == symm_x)
        mirror_x(pixbuf, W, H);
      else
      if (symm == symm_xy)
        mirror_x(pixbuf, W, H - (H >> 1));
      if ((symm == symm_y) || (symm == symm_xy))
        mirror_y(pixbuf, W, H);
      if (symm == symm_point)
        mirror_p(pixbuf, W, H);

      char *png_out_path = NULL;
      if (png_out_dir) {
        static char path[1000];
        snprintf(path, 999, "%s/out%05d.png", png_out_dir, frames_rendered);
        png_out_path = path;
      }
      render(winbuf, winW, winH, &palette, pixbuf, W, H, multiply_pixels, colorshift, out_stream, png_out_path);

      SDL_UpdateTexture(texture, NULL, winbuf, winW * sizeof(Uint32));

      SDL_RenderClear(renderer);
      SDL_RenderCopy(renderer, texture, NULL, NULL);
      SDL_RenderPresent(renderer);

      frames_rendered ++;
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
              underdampen += .0002;
              break;
            case SDLK_LEFT:
              underdampen -= .0002;
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

            case '+':
            case '=':
              if (apex_r < W)
                apex_r ++;
              break;

            case '-':
              if (apex_r > 1)
                apex_r --;
              break;

            case 't':
              underdampen = .996 + .006;
              break;
            case 'r':
              underdampen = .996 + .002;
              break;
            case 'e':
              underdampen = .996;
              break;
            case 'w':
              underdampen = .996 - .002;
              break;
            case 'q':
              underdampen = .996 - .006;
              break;

            default:

              if ((c >= '1') && (c <= '9')) {
                apex_r = 1 + c - '1';
              }
              break;
          }
          }

          printf("underdampen=%f  wavy_amp=%f  symm=%s  apex_r=%d\n",
            underdampen, wavy_amp, symmetry_name[symm], apex_r);
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
  SDL_Quit();
  return 0;
}
