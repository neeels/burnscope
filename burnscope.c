/* burnscope.c
 * (c) 2014 Neels Hofmeyr <neels@hofmeyr.de>
 *
 * This file is part of burnscope, published under the GNU General Public
 * License v3.
 */

#include <math.h>
#include <stdlib.h>
#include <SDL/SDL.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

typedef float pixel_t;

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

void *malloc_check(int len) {
  void *p = malloc(len);
  if (! p) {
    printf("No mem.");
    exit(-1);
  }
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

float rectangle_sum(pixel_t *pixbuf, int W, int H,
                    int x_start, int y_start,
                    int x_end, int y_end,
                    bool wrap_borders) {
  float sum = 0;

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
      sum += *bufpos;
      bufpos ++;
    }
    bufpos += pitch;
  }
  return sum;
}

float surrounding_sum(pixel_t *pixbuf, const int W, const int H,
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
  float sum;
  int pitch = W - rect_w;
  pixel_t *destpos = destbuf + (rect_y * W) + rect_x;
  for (y = rect_y; y < y_end; y++) {
    for (x = rect_x; x < x_end; x++) {
      sum = surrounding_sum(srcbuf, W, H, x, y, apex_r, wrap_borders) / divider;

      *(destpos ++) = sum - truncf(sum);
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
  int x, y;
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
  int x, y;
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
            palette_t *palette, pixel_t *pixbuf, const int W, const int H,
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
  for (y = 0; y < H; y++) {
    pixel_t *pixbuf_y_pos = pixbufpos;
    for (my = 0; my < multiply_pixels; my ++) {
      pixbufpos = pixbuf_y_pos;
      for (x = 0; x < W; x++) {
        int col = (int)((*pixbufpos) * palette->len);
        col += colorshift;
        col %= palette->len;
        Uint32 raw = *((Uint32*) &(palette->colors[col]));
        for (mx = 0; mx < multiply_pixels; mx++) {
          *screenpos = raw;
          screenpos ++;
        }
        pixbufpos ++;
      }
      screenpos += pitch;
    }
  }

  // Unlock if needed
  if (SDL_MUSTLOCK(screen)) 
    SDL_UnlockSurface(screen);

  // Tell SDL to update the whole screen
  SDL_UpdateRect(screen, 0, 0, winW, winH);
}

void seed1(pixel_t *pixbuf, const int W, const int H, int x, int y,
           pixel_t val) {
  x %= W;
  y %= H;
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
  float underdampen = .9978;
  int frame_period = 70;
  bool usage = false;
  bool error = false;
  bool asymmetrical = false;
  bool wrap_borders = true;
  bool start_blank = false;
  symmetry_t symm = symm_x;

  int c;

  while (1) {
    c = getopt(argc, argv, "a:g:m:p:u:AbBh");
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
"  -g WxH  Set animation width and height in number of pixels.\n"
"  -p ms   Set frame period to <ms> milliseconds (slow things down).\n"
"          If zero, run as fast as possible. Default is %d.\n"
"  -m N    Multiply each pixel N times in width and height, to give a larger\n"
"          picture. This will also multiply the window size.\n"
"  -a W    Set apex radius, i.e. the blur distance. Default is %d.\n"
"  -u N.n  Set underdampening factor (decimal). Default is %.3f.\n"
"          Reduces normal blur dampening by this factor.\n"
"  -b      Assume zeros around borders. Default is to wrap around borders.\n"
"  -B      Start out blank. (Use 's' key to plant seeds while running.)\n"
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


  SDL_Surface *screen;

  if ( SDL_Init(SDL_INIT_VIDEO) < 0 ) 
  {
    fprintf(stderr, "Unable to init SDL: %s\n", SDL_GetError());
    exit(1);
  }
  atexit(SDL_Quit);
    
  screen = SDL_SetVideoMode(winW, winH, 32, SDL_SWSURFACE);
  if ( screen == NULL ) 
  {
    fprintf(stderr, "Unable to set %dx%d video: %s\n", winW, winH, SDL_GetError());
    exit(1);
  }
  
  SDL_WM_SetCaption("burnscope", "burnscope");

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
  make_palette(&palette, 4096,
               palette_points, n_palette_points,
               screen->format);

  pixel_t *buf1 = malloc_check(W * H * sizeof(pixel_t));
  pixel_t *buf2 = malloc_check(W * H * sizeof(pixel_t));
  bzero(buf1, W * H * sizeof(pixel_t));
  bzero(buf2, W * H * sizeof(pixel_t));

  pixel_t *pixbuf = buf1;
  pixel_t *swapbuf = buf2;

  int rseed = time(NULL);
  printf("random seed: %d\n", rseed);
  srandom(rseed);

  symmetry_t symm_was = -1;


  if (! start_blank) {
    int i, j;
    j = (W * H / 5) / apex_r;
    for (i = 0; i < j; i ++) {
      seed(pixbuf, W, H, random() % (W), random() % (H), .70, apex_r);
    }
  }

  int last_ticks = SDL_GetTicks() - frame_period;
  bool seed_key_down = false;
  int do_seed = 0;

  float t = 0;
  int cc = 0;
  float wavy = 0;
  int colorshift = 0;

  while (1)
  {
    bool do_render = false;

    t = (float)SDL_GetTicks() / (1200. * 2*M_PI);
    wavy = -.003 * sin(t);
    colorshift = palette.len * (0.5 + 0.5 * cos(t*0.35162748327));
    
    if ((++cc) > 40) {
      printf("%.6f %6d %6d    \r", underdampen + wavy, colorshift, palette.len);
      fflush(stdout);
      cc = 0;
    }

    float divider = 1 + 2 * apex_r;
    divider *= divider * (underdampen + wavy);

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
      while (do_seed) {
        seed(pixbuf, W, H, random() % (W>>1), random() % (H>>1), 200, apex_r);
        do_seed --;
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

      render(screen, winW, winH, &palette, pixbuf, W, H, multiply_pixels, colorshift);
    }
    else
      SDL_Delay(5);

    float underdampen_was = underdampen;
    SDL_Event event;
    while (SDL_PollEvent(&event)) 
    {
      switch (event.type) 
      {
        case SDL_KEYDOWN:
          // If escape is pressed, return (and thus, quit)

          switch(event.key.keysym.sym) {
            case SDLK_ESCAPE:
              return 0;

            case SDLK_RIGHT:
              underdampen += .0002;
              break;
            case SDLK_LEFT:
              underdampen -= .0002;
              break;
            case SDLK_UP:
              underdampen += .000001;
              break;
            case SDLK_DOWN:
              underdampen -= .000001;
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

            default:
              break;
          }
          break;

        case SDL_KEYUP:
          if (event.key.keysym.sym == 's') {
            seed_key_down = false;
          }
          break;

        case SDL_QUIT:
          return 0;
      }
    }

    if (underdampen_was != underdampen) {
      printf("%f\r", underdampen);
      fflush(stdout);
    }

    if (symm != symm_was) {
      printf("%s\n", symmetry_name[symm]);
      symm_was = symm;
    }
  }
  return 0;
}
