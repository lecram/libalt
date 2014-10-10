#ifndef _ALT_H
#define _ALT_H

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

#define ALT_INIT_BULK 7

#define ALT_AT(A, I) (A)->items + (I) * (A)->size

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
    int x0, y0, x1, y1;
    /* Arrays of arrays of alt_cross_t. */
    alt_array_t **hori, **vert, **extr;
} alt_window_t;

void alt_bbmargin(alt_bbox_t *bb, double mg);
void alt_bound(alt_endpt_t *points, int count, bool rounded, alt_bbox_t *bb);

alt_array_t *alt_new_array(size_t item_size, unsigned int init_bulk);
void alt_resize_array(alt_array_t *array, unsigned int min_bulk);
void alt_push(alt_array_t *array, void *item);
void alt_pop(alt_array_t *array, void *item);
void alt_sort(alt_array_t *array, int (*comp)(const void *, const void *));
void alt_del_array(alt_array_t **array);

alt_window_t *alt_new_window(int x0, int y0, int x1, int y1);
void alt_scan(alt_window_t *window, alt_endpt_t *points, int count, double range);
void alt_del_window(alt_window_t **window);

#endif /* _ALT_H */
