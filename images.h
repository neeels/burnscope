
#include "png.h"
#include <sys/types.h>
#include <dirent.h>

typedef struct {
  pixel_t *data;
  int width;
  int height;
  const char *path;
} image_t;


int cmp_str(const void *a, const void *b) {
  return strcmp(*(const char**)a, *(const char **)b);
}


void read_images(const char *images_dir, image_t **images_p, int *n_images_p, int maxW, int maxH) {

  *images_p = NULL;
  *n_images_p = 0;

  int n = 0;
  image_t *images = NULL;
  DIR *dir = opendir(images_dir);


  if (! dir) {
    fprintf(stderr, "Cannot open dir: %s\n", images_dir);
    return;
  }

  // get all file paths and sort alphabetically
  char **files = NULL;
  int n_files = 0;

  while (1) {
    struct dirent *entry = readdir(dir);
    if (! entry)
      break;
    const char *fname = entry->d_name;

    if (strcmp(".", fname) == 0 || strcmp("..", fname) == 0)
      continue;

    int l = strlen(images_dir) + 1 + strlen(fname) + 1;
    char *fpath = malloc(l);
    snprintf(fpath, l, "%s/%s", images_dir, fname);

    if (access(fpath, R_OK) != 0) {
      fprintf(stderr, "cannot read file: %s\n", fpath);
      continue;
    }

    n_files ++;
    files = realloc(files, n_files * sizeof(*files));
    files[n_files - 1] = fpath;
  }

  qsort(files, n_files, sizeof(*files), cmp_str);

  int path_i;
  int images_max_w = 0;
  int images_max_h = 0;

  // get the image sizes first
  for (path_i = 0; path_i < n_files; path_i ++) {
    const char *fpath = files[path_i];

    FILE *infile = fopen(fpath, "r");
    if (infile == NULL) {
      fprintf(stderr, "cannot open file: %s\n", fpath);
      continue;
    }

    uint8_t sig[8];

    fread(sig, 1, 8, infile);
    if (!png_check_sig(sig, 8)) {
      fprintf(stderr, "no png file: %s\n", fpath);
      continue;
    }

    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    png_uint_32 width, height;
    int bit_depth, color_type;


    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) {
      fprintf(stderr, "no mem for png struct while loading file: %s\n", fpath);
      break;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
      png_destroy_read_struct(&png_ptr, NULL, NULL);
      fprintf(stderr, "no mem for png struct while loading file: %s\n", fpath);
      break;
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      fprintf(stderr, "some fucking whatever while loading file: %s\n", fpath);
      continue;
    }

    png_init_io(png_ptr, infile);
    png_set_sig_bytes(png_ptr, 8);

    png_read_info(png_ptr, info_ptr);

    png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type,
      NULL, NULL, NULL);

    printf(" %4dx%4d %s\n", (int)width, (int)height, fpath);
    images_max_w = max(images_max_w, width);
    images_max_h = max(images_max_h, height);

    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
  }

  printf("Largest image dimensions: %dx%d\n", images_max_w, images_max_h);


  int size_shift = 0;
  int want_width = images_max_w;
  int want_height = images_max_h;
  while ((want_width > maxW) || (want_height > maxH)) {
    size_shift ++;
    want_width >>= 1;
    want_height >>= 1;
  }

  if (size_shift)
    printf("Shrinking all images, dividing by %d to fit %dx%d\n",
           1 << size_shift, maxW, maxH);



  for (path_i = 0; path_i < n_files; path_i ++) {
    const char *fpath = files[path_i];

    FILE *infile = fopen(fpath, "r");
    if (infile == NULL) {
      fprintf(stderr, "cannot open file: %s\n", fpath);
      continue;
    }

    uint8_t sig[8];

    fread(sig, 1, 8, infile);
    if (!png_check_sig(sig, 8)) {
      fprintf(stderr, "no png file: %s\n", fpath);
      continue;
    }

    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    png_uint_32 width, height;
    int bit_depth, color_type;
    uint8_t *image_data = NULL;


    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) {
      fprintf(stderr, "no mem for png struct while loading file: %s\n", fpath);
      break;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
      png_destroy_read_struct(&png_ptr, NULL, NULL);
      fprintf(stderr, "no mem for png struct while loading file: %s\n", fpath);
      break;
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      fprintf(stderr, "some fucking whatever while loading file: %s\n", fpath);
      continue;
    }

    png_init_io(png_ptr, infile);
    png_set_sig_bytes(png_ptr, 8);

    png_read_info(png_ptr, info_ptr);

    png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type,
      NULL, NULL, NULL);

    printf(" %4dx%4d", (int)width, (int)height);
    want_width = width >> size_shift;
    want_height = height >> size_shift;
    if (size_shift)
      printf(" --> %4dx%4d", want_width, want_height);
    printf(" %s\n", fpath);


    png_uint_32  i, rowbytes;
    int channels;
    png_bytepp  row_pointers = NULL;

    /* expand palette images to RGB, low-bit-depth grayscale images to 8 bits,
     * transparency chunks to full alpha channel; strip 16-bit-per-sample
     * images to 8 bits per sample; and convert grayscale to RGB[A] */

    if (color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_expand(png_ptr);
    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
        png_set_expand(png_ptr);
    if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
        png_set_expand(png_ptr);
    if (bit_depth == 16)
        png_set_strip_16(png_ptr);
    if (color_type == PNG_COLOR_TYPE_GRAY ||
        color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
        png_set_gray_to_rgb(png_ptr);


    /* unlike the example in the libpng documentation, we have *no* idea where
     * this file may have come from--so if it doesn't have a file gamma, don't
     * do any correction ("do no harm") */

    /*
    if (png_get_gAMA(png_ptr, info_ptr, &gamma))
        png_set_gamma(png_ptr, display_exponent, gamma);
        */


    /* all transformations have been registered; now update info_ptr data,
     * get rowbytes and channels, and allocate image memory */

    png_read_update_info(png_ptr, info_ptr);

    rowbytes = png_get_rowbytes(png_ptr, info_ptr);
    channels = (int)png_get_channels(png_ptr, info_ptr);

    image_data = (uint8_t *)malloc(rowbytes*height);

    if (image_data == NULL) {
        png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
        fprintf(stderr, "no mem for image data while loading file: %s\n", fpath);
        break;
    }
    row_pointers = (png_bytepp)malloc(height*sizeof(png_bytep));
    if (row_pointers == NULL) {
        fprintf(stderr, "no mem for row pointers while loading file: %s\n", fpath);
        png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
        free(image_data);
        image_data = NULL;
        break;
    }


    /* set the individual row_pointers to point at the correct offsets */
    for (i = 0;  i < height;  ++i)
        row_pointers[i] = image_data + i*rowbytes;


    /* now we can go ahead and just read the whole image */
    png_read_image(png_ptr, row_pointers);

    png_read_end(png_ptr, NULL);


    // got the image in mem. convert to internal format.


    n ++;
    images = realloc(images, sizeof(image_t) * n);
    if (! images) {
      fprintf(stderr, "cannot realloc images index to %d\n", (int)(sizeof(image_t) * n));
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      free(row_pointers);
      free(image_data);
      break;
    }

    image_t *img = &images[n - 1];
    img->width = want_width;
    img->height = want_height;
    img->data = malloc(want_width * want_height * sizeof(*(img->data)));
    img->path = fpath;

    if (! img->data) {
      fprintf(stderr, "no mem for image while loading %s\n", fpath);
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      free(row_pointers);
      free(image_data);
      break;
    }

    pixel_t pixmax = 0;
    pixel_t *dest = img->data;

    int src_skip = (1 << size_shift) - 1; // size_shift == 0 ==> src_skip = 0
    int y;
    for (y = 0; y < want_height; y ++) {
      uint8_t *src = row_pointers[y << size_shift];

      int x;
      for (x = 0; x < want_width; x++) {

        int ch;
        pixel_t pixel = 0;

        for (ch = 0; ch < channels; ch++, src++) {
          if (ch == 3) {
            // opacity channel
            pixel *= ((double)(*src)) / 255.;
          }
          else
            pixel += *src;
        }


        *(dest++) = pixel;

        if (pixel > pixmax)
          pixmax = pixel;

        src += src_skip * channels;
      }
    }

    dest = img->data;

    for (y = 0; y < want_height; y ++) {
      int x;
      for (x = 0; x < want_width; x++, dest++) {
        *dest /= pixmax;
      }
    }

    if (row_pointers) {
      free(row_pointers);
      row_pointers = NULL;
    }

    if (image_data) {
        free(image_data);
        image_data = NULL;
    }

    if (png_ptr && info_ptr) {
        png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
        png_ptr = NULL;
        info_ptr = NULL;
    }
  }

  *images_p = images;
  *n_images_p = n;
}
