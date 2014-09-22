#include <SDL/SDL.h>
#include <assert.h>

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

typedef struct {
  unsigned int n_points;
  palette_point_t points[];
} palette_def_t;

palette_def_t pal1 = 
    {
      11,
      {
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
      }
    };

palette_def_t pal2 =
    {
      2,
      {
        { 0, 0, 0, .1 },
        { 0.5 + 3./256, 0, 1, .1 },
      }
    };

#define O .0
#define i (O + .02)
#define g (I * .98)
#define I 1.
#define d -.03

palette_def_t pal3 =
    {
      60,
      {
        {0.+0./6,   O, i, i },
     {d+ 0.+.5/6,   O, i, i },
        {0.+.5/6,   g, I, I },
     {d+ 0.+1./6,   g, I, I },
        {0.+1./6,   O, O, i },
     {d+ 0.+1.5/6,  O, O, i },
        {0.+1.5/6,  g, g, I },
     {d+ 0.+3./6,   g, g, I },
        {0.+3./6,   i, O, O },
     {d+ 0.+3.5/6,  i, O, O },
        {0.+3.5/6,  I, g, g },
     {d+ 0.+4.5/6,  I, g, g },
        {0.+4.5/6,  i, i, O },
     {d+ 0.+4.8/6,  i, i, O },
        {0.+4.8/6,  I, I, g },
     {d+ 0.+5.25/6, I, I, g },
        {0.+5.25/6, O, i, O },
     {d+ 0.+5.7/6,  O, i, O },
        {0.+5.7/6,  g, I, g },
     {d+ 1.+0./6,   g, I, g },
        {1.+0./6,   O, i, i },
     {d+ 1.+.5/6,   O, i, i },
        {1.+.5/6,   g, I, I },
     {d+ 1.+1./6,   g, I, I },
        {1.+1./6,   O, O, i },
     {d+ 1.+1.5/6,  O, O, i },
        {1.+1.5/6,  g, g, I },
     {d+ 1.+3./6,   g, g, I },
        {1.+3./6,   i, O, O },
     {d+ 1.+3.5/6,  i, O, O },
        {1.+3.5/6,  I, g, g },
     {d+ 1.+4.5/6,  I, g, g },
        {1.+4.5/6,  i, i, O },
     {d+ 1.+4.8/6,  i, i, O },
        {1.+4.8/6,  I, I, g },
     {d+ 1.+5.25/6, I, I, g },
        {1.+5.25/6, O, i, O },
     {d+ 1.+5.7/6,  O, i, O },
        {1.+5.7/6,  g, I, g },
     {d+ 2.+0./6,   g, I, g },
        {2.+0./6,   O, i, i },
     {d+ 2.+.5/6,   O, i, i },
        {2.+.5/6,   g, I, I },
     {d+ 2.+1./6,   g, I, I },
        {2.+1./6,   O, O, i },
     {d+ 2.+1.5/6,  O, O, i },
        {2.+1.5/6,  g, g, I },
     {d+ 2.+3./6,   g, g, I },
        {2.+3./6,   i, O, O },
     {d+ 2.+3.5/6,  i, O, O },
        {2.+3.5/6,  I, g, g },
     {d+ 2.+4.5/6,  I, g, g },
        {2.+4.5/6,  i, i, O },
     {d+ 2.+4.8/6,  i, i, O },
        {2.+4.8/6,  I, I, g },
     {d+ 2.+5.25/6, I, I, g },
        {2.+5.25/6, O, i, O },
     {d+ 2.+5.7/6,  O, i, O },
        {2.+5.7/6,  g, I, g },
     {d+ 3,         O, i, i },
      }
    };
#undef O
#undef g
#undef I
#undef i
#undef d

palette_def_t pal4 = 
    {
      4,
      {
        {0, 0, 0, 0},
        {0.05, 1, 1, 1},
        {0.5, 1, 1, 1},
        {0.55, 0, 0, 0},
      }
    };

palette_def_t pal5 =
    {
      2,
      {
        { 0, .3, .3, .3 },
        { 0.5 + 3./256, .3, .3, .3 },
      }
    };


#define n_palettes 5

palette_def_t *palette_defs[n_palettes] = {
    &pal1,
    &pal2,
    &pal3,
    &pal4,
    &pal5,
  };

palette_t palettes[n_palettes];


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
                  palette_def_t *palette_def,
                  SDL_PixelFormat *format) {
  int i;
  int n_points = palette_def->n_points;
  palette_point_t *points = palette_def->points;

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
    float norm_factor = last_p->pos * (n_points + 1) / n_points;
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

void make_palettes(SDL_PixelFormat *format) {
  int pal_i;
  for (pal_i = 0; pal_i < n_palettes; pal_i ++) {
    make_palette(&palettes[pal_i], PALETTE_LEN,
                 palette_defs[pal_i],
                 format);
  }
}

void blend_palettes(palette_t *dst, palette_t *src1, palette_t *src2, float blend) {
  assert((dst->len == src1->len) && (dst->len == src2->len));

  int i, j;
  for (i = 0; i < dst->len; i++) {
    unsigned char *a = (unsigned char*)(&src1->colors[i]);
    unsigned char *b = (unsigned char*)(&src2->colors[i]);
    unsigned char *x = (unsigned char*)(&dst->colors[i]);
    for (j = 0; j < 3; j++) {
      x[j] = (1. - blend) * a[j] + blend * b[j];
    }
  }
}

