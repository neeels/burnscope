
#include "png.h"
#include <sys/types.h>
#include <dirent.h>

typedef struct {
  pixel_t *data;
  int width;
  int height;
} image_t;



void read_images(const char *images_dir, image_t **images_p, int *n_images_p) {

  *images_p = NULL;
  *n_images_p = 0;

  int n = 0;
  image_t *images = NULL;
  DIR *dir = opendir(images_dir);


  if (! dir) {
    fprintf(stderr, "Cannot open dir: %s\n", images_dir);
    return;
  }

  char fpath[PATH_MAX + 1];
  fpath[PATH_MAX] = 0;

  while (1) {
    struct dirent *entry = readdir(dir);
    if (! entry)
      break;
    const char *fname = entry->d_name;
    snprintf(fpath, PATH_MAX, "%s/%s", images_dir, fname);

    if (strcmp(".", fname) == 0 || strcmp("..", fname) == 0)
      continue;

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

    printf("%s %dx%d\n", fname, (int)width, (int)height);

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
    img->width = width;
    img->height = height;
    img->data = malloc(width * height * sizeof(*(img->data)));

    if (! img->data) {
      fprintf(stderr, "no mem for image while loading %s\n", fpath);
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      free(row_pointers);
      free(image_data);
      break;
    }

    pixel_t pixmax = 0;
    pixel_t *dest = img->data;

    int row_i;
    for (row_i = 0; row_i < height; row_i ++) {
      uint8_t *src = row_pointers[row_i];

      int x;
      for (x = 0; x < width; x++) {

        int ch;
        pixel_t pixel = 0;

        for (ch = 0; ch < channels; ch++, src++) {
          pixel += *src;
          if (ch == 3) {
            // opacity channel
            pixel *= ((double)(*src)) / 255.;
          }
        }


        *(dest++) = pixel;
        if (pixel > pixmax)
          pixmax = pixel;
      }
    }

    dest = img->data;

    for (row_i = 0; row_i < height; row_i ++) {
      int x;
      for (x = 0; x < width; x++, dest++) {
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
