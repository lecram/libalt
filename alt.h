#ifndef _ALT_H
#define _ALT_H

#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>

#define ALT_MIN(A, B) ((A) < (B) ? (A) : (B))
#define ALT_MAX(A, B) ((A) > (B) ? (A) : (B))

/* Bounding box. */
typedef struct {
    double x0, y0, x1, y1;
} alt_bbox_t;

/* End point for line segments. */
typedef struct {
    double x, y;
} alt_endpt_t;

/* Control point for Bézier curves. */
typedef struct {
    double x, y;
    int on;
} alt_ctrpt_t;

#define ALT_PIX(I, X, Y) ((I)->data[4*((Y)*(I)->width+(X))])

/* Raw 32-bit RGBA image. */
typedef struct {
    int width, height;  /* image size in pixels */
    uint8_t *data;     /* width*height pixels in format 0xRRGGBBAA */
} alt_image_t;

#define ALT_INIT_BULK 7

#define ALT_AT(A, I) ((A)->items + (I) * (A)->size)
#define ALT_CAT(T, A, I) ((T *) ALT_AT((A), (I)))

/* Generic dynamic array. */
typedef struct {
    unsigned int bulk;  /* current capacity */
    unsigned int count; /* current ocupation */
    size_t size;        /* item size in bytes */
    void *items;        /* actual array contents */
} alt_array_t;

/* Signed scanline intersection. */
typedef struct {
    double dist;    /* distance from scanline origin */
    int sign;       /* 1 for off-on crossing, -1 for on-off */
} alt_cross_t;

/* Scan window. */
typedef struct {
    int x0, y0, width, height;
    /* Arrays of arrays of alt_cross_t. */
    alt_array_t **hori, **vert, **extr;
} alt_window_t;

bool alt_bbisnull(alt_bbox_t *bb);
void alt_bbsetnull(alt_bbox_t *bb);
void alt_bbinter(alt_bbox_t *bb1, alt_bbox_t *bb2);
void alt_bbunion(alt_bbox_t *bb1, alt_bbox_t *bb2);
void alt_bbmargin(alt_bbox_t *bb, double mg);
void alt_bound(alt_endpt_t *points, int count, alt_bbox_t *bb);

uint32_t alt_pack_color(uint8_t r, uint8_t g, uint8_t b, uint8_t a);
void alt_unpack_color(uint32_t color, uint8_t *r, uint8_t *g, uint8_t *b, uint8_t *a);
alt_image_t *alt_new_image(int width, int height);
void alt_get_pixel(alt_image_t *image, int x, int y,
                   uint8_t *r, uint8_t *g, uint8_t *b, uint8_t *a);
void alt_set_pixel(alt_image_t *image, int x, int y,
                   uint8_t r, uint8_t g, uint8_t b, uint8_t a);
void alt_clear(alt_image_t *image, uint32_t color);
void alt_blend(alt_image_t *image, int x, int y,
               uint8_t r, uint8_t g, uint8_t b, uint8_t a);
void alt_save_pam(alt_image_t *image, char *fname);
void alt_del_image(alt_image_t **image);

alt_array_t *alt_new_array(size_t item_size, unsigned int init_bulk);
void alt_resize_array(alt_array_t *array, unsigned int min_bulk);
void alt_push(alt_array_t *array, void *item);
void alt_pop(alt_array_t *array, void *item);
void alt_sort(alt_array_t *array, int (*comp)(const void *, const void *));
void alt_del_array(alt_array_t **array);

int alt_comp_cross(const void *a, const void *b);

alt_window_t *alt_new_window(alt_bbox_t *bb);
void alt_scan(alt_window_t *window, alt_endpt_t *points, int count, double range);
void alt_windredux(alt_window_t *window);
void alt_del_window(alt_window_t **window);

double alt_scanrange(alt_array_t *scanline, double x);
double alt_dist(alt_window_t *window, double x, double y, double r);

void alt_draw(alt_image_t *image, alt_window_t *window,
              uint32_t fill, uint32_t strk, double thick);

#endif /* _ALT_H */
