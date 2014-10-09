#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "alt.h"

/* Add a margin to the bounding box at `bb`. */
void
alt_bbmargin(alt_bbox_t *bb, double mg)
{
    bb->x0 -= mg;
    bb->y0 -= mg;
    bb->x1 += mg;
    bb->y1 += mg;
}

/* Compute the bounding box of `count` points at `points`.
 * If `rounded` is nonzero, the result is fitted to the integer lattice.
 * The result is written to `bb`.
 */
void
alt_bound(alt_endpt_t *points, int count, bool rounded, alt_bbox_t *bb)
{
    int i;
    bb->x0 = bb->x1 = points[0].x;
    bb->y0 = bb->y1 = points[0].y;
    for (i = 1; i < count; i++) {
        if (points[i].x < bb->x0)
            bb->x0 = points[i].x;
        else if (points[i].x > bb->x1)
            bb->x1 = points[i].x;
        if (points[i].y < bb->y0)
            bb->y0 = points[i].y;
        else if (points[i].y > bb->y1)
            bb->y1 = points[i].y;
    }
    if (rounded) {
        bb->x0 = floor(bb->x0);
        bb->y0 = floor(bb->y0);
        bb->x1 = ceil(bb->x1);
        bb->y1 = ceil(bb->y1);
    }
}

/* Create a new array. Return NULL if there is not enough memory. */
alt_array_t *
alt_new_array(size_t item_size, unsigned int init_bulk)
{
    alt_array_t *array;
    array = (alt_array_t *) malloc(sizeof(alt_array_t));
    if (array == NULL) return NULL;
    array->bulk = init_bulk;
    array->count = 0;
    array->size = item_size;
    array->items = malloc(array->bulk * array->size);
    if (array->items == NULL) {
        free(array);
        return NULL;
    }
    return array;
}

/* Resize the array to accomodate at least `min_bulk` items. */
void
alt_resize_array(alt_array_t *array, unsigned int min_bulk)
{
    while (array->bulk > min_bulk)
        array->bulk >>= 1;
    while (array->bulk < min_bulk)
        array->bulk <<= 1;
    array->items = realloc(array->items, array->bulk * array->size);
}

/* Add item to the end of `array`. */
void
alt_push(alt_array_t *array, void *item)
{
    if (array->count == array->bulk)
        alt_resize_array(array, array->count + 1);
    memcpy(ALT_AT(array, array->count), item, array->size);
    array->count++;
}

/* Remove item from the end of `array`.
 * The item removed is copied into `item`, if it's not null.
 */
void
alt_pop(alt_array_t *array, void *item)
{
    array->count--;
    if (item != NULL)
        memcpy(item, ALT_AT(array, array->count), array->size);
    if (array->count < array->bulk >> 1)
        alt_resize_array(array, array->count);
}

/* Sort array according to `comp` (see stdlib's qsort()). */
void
alt_sort(alt_array_t *array, int (*comp)(const void *, const void *))
{
    qsort(array->items, (size_t) array->count, array->size, comp);
}

/* Delete array and its content from memory. */
void
alt_del_array(alt_array_t **array)
{
    free((*array)->items);
    free(*array);
    *array = NULL;
}

/* Create a new scan window. */
alt_window_t *
alt_new_window(int x0, int y0, int x1, int y1)
{
    alt_window_t *window;
    int width, height;
    window = (alt_window_t *) malloc(sizeof(alt_window_t));
    if (window == NULL) return NULL;
    window->x0 = x0; window->y0 = y0;
    window->x1 = x1; window->y1 = y1;
    width  = x1-x0+1;
    height = y1-y0+1;
    window->hori = alt_new_array(sizeof(alt_cross_t), height);
    window->vert = alt_new_array(sizeof(alt_cross_t), width);
    window->extr = alt_new_array(sizeof(alt_cross_t), height);
    return window;
}
