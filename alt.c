#include <stdlib.h>
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
    array->bulk = init_bulk || ALT_INIT_BULK;
    array->count = 0;
    array->size = item_size;
    array->items = malloc(array->bulk * array->size);
    if (array->items == NULL) return NULL;
    return array;
}

/* Resize the array to accomodate at least `min_bulk` items. */
void
alt_resize_array(alt_array_t *array, unsigned int min_bulk)
{
    min_bulk = min_bulk || 2;
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

/* Delete `array` and its content from memory. */
void
alt_del_array(alt_array_t **array)
{
    free((*array)->items);
    free(*array);
    *array = NULL;
}

/* Intersection comparison for qsort()-like interfaces. */
int
alt_comp_cross(const void *a, const void *b)
{
    const alt_cross_t *arg1 = a;
    const alt_cross_t *arg2 = b;
    return arg1->dist - arg2->dist;
}

/* Create a new scan window. */
alt_window_t *
alt_new_window(alt_bbox_t *bb)
{
    alt_window_t *window;
    int i;
    window = (alt_window_t *) malloc(sizeof(alt_window_t));
    if (window == NULL) return NULL;
    window->x0 = (int) bb->x0;
    window->y0 = (int) bb->y0;
    window->width  = (int) (bb->x1 - bb->x0 + 1);
    window->height = (int) (bb->y1 - bb->y0 + 1);
    window->hori = (alt_array_t **) malloc(window->height * sizeof(alt_array_t *));
    window->vert = (alt_array_t **) malloc(window->width  * sizeof(alt_array_t *));
    window->extr = (alt_array_t **) malloc(window->height * sizeof(alt_array_t *));
    if (window->hori == NULL || window->vert == NULL || window->extr == NULL)
        return NULL;
    for (i = 0; i < window->width; i++) {
        window->vert[i] = alt_new_array(sizeof(alt_cross_t), ALT_INIT_BULK);
        if (window->vert[i] == NULL) return NULL;
    }
    for (i = 0; i < window->height; i++) {
        window->hori[i] = alt_new_array(sizeof(alt_cross_t), ALT_INIT_BULK);
        window->extr[i] = alt_new_array(sizeof(alt_cross_t), ALT_INIT_BULK);
        if (window->hori[i] == NULL || window->extr[i] == NULL) return NULL;
    }
    return window;
}

/* Scan `count` points at `points` and add intersections to `window`.
 * `range` is the line width for the extra scan.
 */
void
alt_scan(alt_window_t *window, alt_endpt_t *points, int count, double range)
{
    alt_endpt_t pa, pb;
    alt_array_t *scanline;
    alt_cross_t cross;
    double radius;
    double hxa, hya, hxb, hyb;
    double vxa, vya, vxb, vyb;
    double ex0, ey0, ex1, ey1;
    double slope;
    int hsign, vsign;
    int sx, sy;
    int i, y;
    radius = range / 2;
    for (i = 0; i < count-1; i++) {
        pa = points[i];
        pb = points[i+1];
        /*  Define a segment from (vxa, vya) to (vxb, vyb) that coincides with 
         * pa-pb, but is guaranteed to be x-ascending.
         */
        if (pa.x < pb.x) {
            hsign = +1;
            vxa = pa.x; vya = pa.y; vxb = pb.x; vyb = pb.y;
        }
        else {
            hsign = -1;
            vxa = pb.x; vya = pb.y; vxb = pa.x; vyb = pa.y;
        }
        /*  Define a segment from (hxa, hya) to (hxb, hyb) that coincides with 
         * pa-pb, but is guaranteed to be y-ascending.
         */
        if (pa.y < pb.y) {
            vsign = +1;
            hxa = pa.x; hya = pa.y; hxb = pb.x; hyb = pb.y;
        }
        else {
            vsign = -1;
            hxa = pb.x; hya = pb.y; hxb = pa.x; hyb = pa.y;
        }
        vxa = round(vxa); vxb = round(vxb);
        hya = round(hya); hyb = round(hyb);
        /*  Add orthogonal intersections to extra scan. This is suboptimal for
         * diagonal segments, specially when fabs(slope) ~ 1.
         */
        ex0 = floor(vxa-radius); ey0 = floor(hya-radius);
        ex1 = ceil(vxb+radius)+1; ey1 = ceil(hyb+radius)+1;
        cross.dist = ex0; cross.sign = +1;
        for (y = (int) ey0; y < (int) ey1; y++) {
            scanline = window->extr[y-window->y0];
            alt_push(scanline, &cross);
        }
        cross.dist = ex1; cross.sign = -1;
        for (y = (int) ey0; y < (int) ey1; y++) {
            scanline = window->extr[y-window->y0];
            alt_push(scanline, &cross);
        }
        if (pa.y == pb.y) {
            /* Horizontal segment. */
            cross.dist = pa.y; cross.sign = hsign;
            while (vxa < vxb) {
                scanline = window->vert[((int) vxa)-window->x0];
                alt_push(scanline, &cross);
                vxa++;
            }
        }
        else if (pa.x == pb.x) {
            /* Vertical segment. */
            cross.dist = pa.x; cross.sign = vsign;
            while (hya < hyb) {
                scanline = window->hori[((int) hya)-window->y0];
                alt_push(scanline, &cross);
                hya++;
            }
        }
        else {
            /* Diagonal segment. */
            sy = vya < vyb ? +1 : -1;
            slope = (vxb-vxa) / (hyb-hya);
            cross.sign = hsign;
            while (vxa < vxb) {
                scanline = window->vert[((int) vxa)-window->x0];
                cross.dist = vya;
                alt_push(scanline, &cross);
                vya += sy / slope;
                vxa++;
            }
            sx = hxa < hxb ? +1 : -1;
            slope = 1 / slope;
            cross.sign = vsign;
            while (hya < hyb) {
                scanline = window->hori[((int) hya)-window->y0];
                cross.dist = hxa;
                alt_push(scanline, &cross);
                hxa += sx / slope;
                hya++;
            }
        }
    }
}

/* Sort intersections and reduce them according to the running winding number.
 * Only crossings that goes to or come from a zero winding number are kept.
 */
void
alt_windredux(alt_window_t *window)
{
    alt_array_t **scans[3], *scanline;
    alt_cross_t cross, *pcross;
    int count[3];
    int winda, windb;
    int i, j, k;
    scans[0] = window->vert;
    scans[1] = window->hori;
    scans[2] = window->extr;
    count[0] = window->width;
    count[1] = count[2] = window->height;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < count[i]; j++) {
            alt_sort(scans[i][j], alt_comp_cross);
            scanline = alt_new_array(sizeof(alt_cross_t), scans[i][j]->count);
            cross.dist = -HUGE_VAL;
            alt_push(scanline, &cross);
            winda = windb = 0;
            for (k = 0; k < (int) scans[i][j]->count; k++) {
                pcross = ALT_CAT(alt_cross_t, scans[i][j], k);
                windb += pcross->sign;
                if (!winda || !windb) {
                    if (pcross->dist == cross.dist) {
                        /* Duplicate crossings cancel each other. */
                        alt_pop(scanline, NULL);
                    }
                    else {
                        cross.dist = pcross->dist;
                        alt_push(scanline, &cross);
                    }
                }
                winda = windb;
            }
            cross.dist = HUGE_VAL;
            alt_push(scanline, &cross);
            alt_resize_array(scanline, scanline->count);
            alt_del_array(&scans[i][j]);
            scans[i][j] = scanline;
        }
    }
}

/* Delete `window` and its content from memory. */
void
alt_del_window(alt_window_t **window)
{
    int i;
    for (i = 0; i < (*window)->width; i++) {
        alt_del_array(&(*window)->vert[i]);
    }
    for (i = 0; i < (*window)->height; i++) {
        alt_del_array(&(*window)->hori[i]);
        alt_del_array(&(*window)->extr[i]);
    }
    free((*window)->hori);
    free((*window)->vert);
    free((*window)->extr);
    free(*window);
    *window = NULL;
}

/* Return the minimum distance from `x` to its closest point on `scanline`.
 * Note that x here is a 1D coordinate that follows the scanline's direction.
 * This function is only used as a helper to alt_dist().
 */
double
alt_scanrange(alt_array_t *scanline, double x)
{
    alt_cross_t *pcross;
    pcross = (alt_cross_t *) scanline->items;
    while ((++pcross)->dist <= x);
    return ALT_MIN(pcross->dist - x, x - (pcross-1)->dist);
}

/* Return the minimum distance between (x, y) and the path scanned into
 * `window`. If that path doesn't intersect the disk of radius `r`
 * centered at (x, y), then return `r`.
 * This function implements the L1L2 algorithm.
 */
double
alt_dist(alt_window_t *window, double x, double y, double r)
{
    double mind, mag;
    double j, d, d1, d2;
    int i, xi, yi;
    xi = x - window->x0 + 1;
    yi = y - window->y0 + 1;
    i = 0;
    mind = r*r;
    mag = r/sqrt(2.0);
    do {
        d1 = i*i;
        j = HUGE_VAL;
        if (xi-i >= 0 && xi-i < window->width)
            j = ALT_MIN(j, alt_scanrange(window->vert[xi-i], y));
        if (xi+i >= 0 && xi+i < window->width)
            j = ALT_MIN(j, alt_scanrange(window->vert[xi+i], y));
        if (yi-i >= 0 && yi-i < window->height)
            j = ALT_MIN(j, alt_scanrange(window->hori[yi-i], x));
        if (yi+i >= 0 && yi+i < window->height)
            j = ALT_MIN(j, alt_scanrange(window->hori[yi+i], x));
        d2 = j*j;
        d = d1 + d2;
        if (d < mind) {
            mag = sqrt(d/2);
            mind = d;
        }
        i++;
    } while (i > (int) mag);
    return mag * sqrt(2);
}
