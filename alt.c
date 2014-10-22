#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "alt.h"

bool
alt_bbisnull(alt_bbox_t *bb)
{
    return bb->x0 >= bb->x1 || bb->y0 >= bb->y1;
}

void
alt_bbsetnull(alt_bbox_t *bb)
{
    bb->x0 = bb->y0 = 0;
    bb->x1 = bb->y1 = -1;
}

/* Compute the intersection of `bb1` and `bb2` and write it to `bb1`. */
void
alt_bbinter(alt_bbox_t *bb1, alt_bbox_t *bb2)
{
    if (alt_bbisnull(bb1)) return;
    if (alt_bbisnull(bb2)) {
        alt_bbsetnull(bb1);
    }
    else {
        bb1->x0 = ALT_MAX(bb1->x0, bb2->x0);
        bb1->y0 = ALT_MAX(bb1->y0, bb2->y0);
        bb1->x1 = ALT_MIN(bb1->x1, bb2->x1);
        bb1->y1 = ALT_MIN(bb1->y1, bb2->y1);
    }
}

/* Compute the union of `bb1` and `bb2` and write it to `bb1`. */
void
alt_bbunion(alt_bbox_t *bb1, alt_bbox_t *bb2)
{
    if (alt_bbisnull(bb1)) {
        bb1->x0 = bb2->x0;
        bb1->y0 = bb2->y0;
        bb1->x1 = bb2->x1;
        bb1->y1 = bb2->y1;
    }
    else if (!alt_bbisnull(bb2)) {
        bb1->x0 = ALT_MIN(bb1->x0, bb2->x0);
        bb1->y0 = ALT_MIN(bb1->y0, bb2->y0);
        bb1->x1 = ALT_MAX(bb1->x1, bb2->x1);
        bb1->y1 = ALT_MAX(bb1->y1, bb2->y1);
    }
}

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
 * The result is written to `bb`.
 */
void
alt_bound(alt_endpt_t *points, int count, alt_bbox_t *bb)
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
}

/* Pack color components into 32-bit value. */
uint32_t
alt_pack_color(uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
    return (r << 24) + (g << 16) + (b << 8) + a;
}

/* Unpack color components from 32-bit value. */
void
alt_unpack_color(uint32_t color, uint8_t *r, uint8_t *g, uint8_t *b, uint8_t *a)
{
    *r = color >> 24;
    *g = (color >> 16) & 0xFF;
    *b = (color >> 8) & 0xFF;
    *a = color & 0xFF;
}

/* Create a new image. Return NULL if there is not enough memory. */
alt_image_t *
alt_new_image(int width, int height)
{
    alt_image_t *image;
    image = (alt_image_t *) malloc(sizeof(alt_image_t));
    if (image == NULL) return NULL;
    image->width  = width;
    image->height = height;
    image->data = (uint8_t *) calloc(4*width*height, sizeof(uint8_t));
    if (image->data == NULL) return NULL;
    return image;
}

/* Get pixel color at (x, y). */
void
alt_get_pixel(alt_image_t *image, int x, int y,
              uint8_t *r, uint8_t *g, uint8_t *b, uint8_t *a)
{
    int i;
    i = 4*(y*image->width+x);
    *r = image->data[i];
    *g = image->data[i+1];
    *b = image->data[i+2];
    *a = image->data[i+3];
}

/* Set pixel color at (x, y). */
void
alt_set_pixel(alt_image_t *image, int x, int y,
              uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
    int i;
    i = 4*(y*image->width+x);
    image->data[i]   = r;
    image->data[i+1] = g;
    image->data[i+2] = b;
    image->data[i+3] = a;
}

/* Fill the entire image with `color`. */
void
alt_clear(alt_image_t *image, uint32_t color)
{
    int i;
    uint8_t r, g, b, a;
    alt_unpack_color(color, &r, &g, &b, &a);
    for (i = 0; i < 4*image->width*image->height; i += 4) {
        image->data[i]   = r;
        image->data[i+1] = g;
        image->data[i+2] = b;
        image->data[i+3] = a;
    }
}

/* Alpha-blend pixel at (x, y) with (r, g, b, a). */
void
alt_blend(alt_image_t *image, int x, int y,
          uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
    uint8_t r0, g0, b0, a0; /* background */
    uint8_t r2, g2, b2, a2; /* blended */
    double da0, da1, da2, o;
    alt_get_pixel(image, x, y, &r0, &g0, &b0, &a0);
    da0 = ((double) a0) / 255;
    da1 = ((double) a) / 255;
    o = da0*(1-da1);
    da2 = da1 + o;
    if (da2 == 0)
        alt_set_pixel(image, x, y, 0, 0, 0, 0);
    else {
        r2 = lround((r*da1 + r0*o) / da2);
        g2 = lround((g*da1 + g0*o) / da2);
        b2 = lround((b*da1 + b0*o) / da2);
        a2 = lround(da2*255);
        alt_set_pixel(image, x, y, r2, g2, b2, a2);
    }
}

/* Save image as PAM file. */
void
alt_save_pam(alt_image_t *image, char *fname)
{
    FILE *pam = fopen(fname, "w");
    if (pam == NULL) return;
    fprintf(pam, "P7\n");
    fprintf(pam, "WIDTH %d\n", image->width);
    fprintf(pam, "HEIGHT %d\n", image->height);
    fprintf(pam, "DEPTH 4\n");
    fprintf(pam, "MAXVAL 255\n");
    fprintf(pam, "TUPLTYPE RGB_ALPHA\n");
    fprintf(pam, "ENDHDR\n");
    fwrite(image->data, sizeof(uint32_t), image->width*image->height, pam);
    fclose(pam);
}

/* Delete `image` and its content from memory. */
void
alt_del_image(alt_image_t **image)
{
    free((*image)->data);
    free(*image);
    *image = NULL;
}

/* Create a new array. Return NULL if there is not enough memory. */
alt_array_t *
alt_new_array(size_t item_size, unsigned int init_bulk)
{
    alt_array_t *array;
    array = (alt_array_t *) malloc(sizeof(alt_array_t));
    if (array == NULL) return NULL;
    array->bulk = init_bulk ? init_bulk : ALT_INIT_BULK;
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
    while (array->bulk > min_bulk)
        array->bulk >>= 1;
    while (array->bulk < min_bulk)
        array->bulk <<= 1;
    array->bulk = array->bulk ? array->bulk : 1;
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
    window->x0 = (int) floor(bb->x0);
    window->y0 = (int) floor(bb->y0);
    window->width  = (int) (ceil(bb->x1) - floor(bb->x0) + 1);
    window->height = (int) (ceil(bb->y1) - floor(bb->y0) + 1);
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

/* Scan segment pa-pb and add intersections to `window`.
 * `range` is the line width for the extra scan.
 */
void
alt_scan(alt_window_t *window, alt_endpt_t *pa, alt_endpt_t *pb, double range)
{
    alt_array_t *scanline;
    alt_cross_t cross;
    double radius;
    double hxa, hya, hxb, hyb;
    double vxa, vya, vxb, vyb;
    double ex0, ey0, ex1, ey1;
    double slope;
    int hsign, vsign;
    int sx, sy;
    int y;
    radius = range / 2;
    /*  Define a segment from (vxa, vya) to (vxb, vyb) that coincides with 
     * pa-pb, but is guaranteed to be x-ascending.
     */
    if (pa->x < pb->x) {
        hsign = +1;
        vxa = pa->x; vya = pa->y; vxb = pb->x; vyb = pb->y;
    }
    else {
        hsign = -1;
        vxa = pb->x; vya = pb->y; vxb = pa->x; vyb = pa->y;
    }
    /*  Define a segment from (hxa, hya) to (hxb, hyb) that coincides with 
     * pa-pb, but is guaranteed to be y-ascending.
     */
    if (pa->y < pb->y) {
        vsign = +1;
        hxa = pa->x; hya = pa->y; hxb = pb->x; hyb = pb->y;
    }
    else {
        vsign = -1;
        hxa = pb->x; hya = pb->y; hxb = pa->x; hyb = pa->y;
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
    if (pa->y == pb->y) {
        /* Horizontal segment. */
        cross.dist = pa->y; cross.sign = hsign;
        while (vxa < vxb) {
            scanline = window->vert[((int) vxa)-window->x0];
            alt_push(scanline, &cross);
            vxa++;
        }
    }
    else if (pa->x == pb->x) {
        /* Vertical segment. */
        cross.dist = pa->x; cross.sign = vsign;
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

/* Scan `count` points at `points` and add intersections to `window`.
 * `range` is the line width for the extra scan.
 */
void
alt_scan_array(alt_window_t *window, alt_endpt_t *points, int count, double range)
{
    int i;
    alt_endpt_t pa, pb;
    for (i = 0; i < count-1; i++) {
        pa = points[i];
        pb = points[i+1];
        alt_scan(window, &pa, &pb, range);
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
    xi = x - window->x0;
    yi = y - window->y0;
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
    } while (i <= (int) mag);
    return mag * sqrt(2);
}

/* Draw scanned figure to `image`. 
 * `fill` is the fill color; `strk` is the stroke color.
 * `thick` is the line width (thickness) for stroke.
 */
void
alt_draw(alt_image_t *image, alt_window_t *window,
         uint32_t fill, uint32_t strk, double thick)
{
    alt_array_t *esl, *hsl;
    alt_cross_t *ecross, *hcross;
    bool border, inside;
    double hlwp, d, m, aa;
    uint8_t fr, fg, fb, fa;
    uint8_t sr, sg, sb, sa;
    int x0, y0, x1, y1, x, y, i;
    alt_unpack_color(fill, &fr, &fg, &fb, &fa);
    alt_unpack_color(strk, &sr, &sg, &sb, &sa);
    hlwp = thick/2 + 0.5;
    x0 = ALT_MAX(0, window->x0);
    y0 = ALT_MAX(0, window->y0);
    x1 = ALT_MIN(image->width, window->x0 + window->width);
    y1 = ALT_MIN(image->height, window->y0 + window->height);
    i = y0 - window->y0;
    for (y = y0; y < y1; y++, i++) {
        esl = window->extr[i];
        hsl = window->hori[i];
        border = inside = false;
        ecross = (alt_cross_t *) esl->items;
        hcross = (alt_cross_t *) hsl->items;
        ecross++; hcross++;
        for (x = window->x0; x < window->x0 + window->width; x++) {
            if (x >= ecross->dist) {
                border = !border;
                ecross++;
            }
            if (x >= hcross->dist) {
                inside = !inside;
                hcross++;
            }
            if (x0 <= x && x < x1 && y0 <= y && y < y1) {
                if (border) {
                    d = alt_dist(window, x, y, hlwp);
                    m = fabs(d);
                    if (inside) {
                        aa = ALT_MIN(d, 1);
                        alt_blend(image, x, y, fr, fg, fb, lround(aa*fa));
                    }
                    if (m < hlwp) {
                        aa = hlwp - m;
                        aa = ALT_MIN(aa, 1);
                        alt_blend(image, x, y, sr, sg, sb, lround(aa*sa));
                    }
                }
                else if (inside) {
                    alt_blend(image, x, y, fr, fg, fb, fa);
                }
            }
        }
    }
}

/* Convert `curve` to polyline and add resulting points to array.
 * `endpts` is an array of alt_endpt_t.
 */
void
alt_add_curve(alt_array_t *endpts, alt_curve_t *curve)
{
    /*  We split a Bézier curve into two subcurves repeatedly until we reach
     * almost-straight curves that are then output as straight segments.
     *  Since the segments must be added in the correct order, we need to do
     * a depth-first search as we create a tree of subcurves.
     */
    alt_array_t *stack;
    alt_curve_t subcurve;
    alt_endpt_t a, b, c;
    alt_endpt_t d, e, f;
    double h;
    stack = alt_new_array(sizeof(alt_curve_t), ALT_INIT_BULK);
    alt_push(stack, curve);
    while (stack->count) {
        alt_pop(stack, &subcurve);
        a = subcurve.a;
        b = subcurve.b;
        c = subcurve.c;
        h = fabs((a.x-c.x)*(b.y-a.y)-(a.x-b.x)*(c.y-a.y)) /
            hypot(c.x-a.x, c.y-a.y);
        if (h > 1) {
            /* Split curve. */
            d.x = (a.x+b.x)/2; d.y= (a.y+b.y)/2;
            f.x = (b.x+c.x)/2; f.y= (b.y+c.y)/2;
            e.x = (d.x+f.x)/2; e.y= (d.y+f.y)/2;
            subcurve.a = e;
            subcurve.b = f;
            subcurve.c = c;
            alt_push(stack, &subcurve);
            subcurve.a = a;
            subcurve.b = d;
            subcurve.c = e;
            alt_push(stack, &subcurve);
        }
        else {
            /* Add point to polyline. */
            alt_push(endpts, &c);
        }
    }
    alt_del_array(&stack);
}

/* Convert a sequence of linked Bézier curves to a polyline.
 *  `ctrpts` is an array of alt_ctrpt_t. The boolean field `.on` indicates
 * whether a point is on-curve or off-curve. On-curve points are end points;
 * off-curve points are Bézier control points. If two off-curve points appear
 * next to each other, an on-curve point is implied halfway between them.
 * Return an array of alt_endpt_t. */
alt_array_t *
alt_unfold(alt_array_t *ctrpts)
{
    alt_ctrpt_t *s, *a, *b, m;
    alt_array_t *endpts;
    alt_curve_t curve;
    int i;
    m.on = true;
    endpts = alt_new_array(sizeof(alt_endpt_t), ALT_INIT_BULK);
    s = ALT_CAT(alt_ctrpt_t, ctrpts, 0);
    b = ALT_CAT(alt_ctrpt_t, ctrpts, 1);
    alt_push(endpts, s);
    if (b->on) {
        alt_push(endpts, b);
        s = b;
    }
    for (i = 2; i < (int) ctrpts->count; i++) {
        a = b;
        b = ALT_CAT(alt_ctrpt_t, ctrpts, i);
        if (a->on) {
            if (b->on) {
                alt_push(endpts, b);
                s = b;
            }
        }
        else {
            if (b->on) {
                curve.a.x = s->x; curve.a.y = s->y;
                curve.b.x = a->x; curve.b.y = a->y;
                curve.c.x = b->x; curve.c.y = b->y;
                alt_add_curve(endpts, &curve);
                s = b;
            }
            else {
                curve.a.x = s->x; curve.a.y = s->y;
                curve.b.x = a->x; curve.b.y = a->y;
                m.x = (a->x + b->x) / 2; m.y = (a->y + b->y) / 2;
                curve.c.x = m.x; curve.c.y = m.y;
                alt_add_curve(endpts, &curve);
                s = &m;
            }
        }
    }
    return endpts;
}

/* Create a bezigon that approximates a circle of center (x, y) and radius r. */
alt_array_t *
alt_circle(double x, double y, double r)
{
    double a;
    alt_ctrpt_t ctrpt;
    alt_array_t *ctrpts;
    a = r * (sqrt(2) - 1);
    ctrpts = alt_new_array(sizeof(alt_ctrpt_t), 10);
    ctrpt.x = x; ctrpt.y = y+r; ctrpt.on = true;
    alt_push(ctrpts, &ctrpt);
    ctrpt.x = x+a; ctrpt.y = y+r; ctrpt.on = false;
    alt_push(ctrpts, &ctrpt);
    ctrpt.x = x+r; ctrpt.y = y+a; ctrpt.on = false;
    alt_push(ctrpts, &ctrpt);
    ctrpt.x = x+r; ctrpt.y = y-a; ctrpt.on = false;
    alt_push(ctrpts, &ctrpt);
    ctrpt.x = x+a; ctrpt.y = y-r; ctrpt.on = false;
    alt_push(ctrpts, &ctrpt);
    ctrpt.x = x-a; ctrpt.y = y-r; ctrpt.on = false;
    alt_push(ctrpts, &ctrpt);
    ctrpt.x = x-r; ctrpt.y = y-a; ctrpt.on = false;
    alt_push(ctrpts, &ctrpt);
    ctrpt.x = x-r; ctrpt.y = y+a; ctrpt.on = false;
    alt_push(ctrpts, &ctrpt);
    ctrpt.x = x-a; ctrpt.y = y+r; ctrpt.on = false;
    alt_push(ctrpts, &ctrpt);
    ctrpt.x = x; ctrpt.y = y+r; ctrpt.on = true;
    alt_push(ctrpts, &ctrpt);
    return ctrpts;
}

/* Reset matrix to identity map. */
void
alt_reset(alt_matrix_t *mat)
{
    mat->a = 1; mat->c = 0; mat->e = 0;
    mat->b = 0; mat->d = 1; mat->f = 0;
}

/* Add a custom transformation to matrix. */
void
alt_add_custom(alt_matrix_t *mat, double a, double b,
               double c, double d, double e, double f)
{
    double na, nb, nc, nd, ne, nf;
    na = a*mat->a + c*mat->b;
    nb = b*mat->a + d*mat->b;
    nc = a*mat->c + c*mat->d;
    nd = b*mat->c + d*mat->d;
    ne = a*mat->e + c*mat->f + e;
    nf = b*mat->e + d*mat->f + f;
    mat->a = na; mat->c = nc; mat->e = ne;
    mat->b = nb; mat->d = nd; mat->f = nf;
}

/* Add a squeeze transformation to matrix. */
void
alt_add_squeeze(alt_matrix_t *mat, double k)
{
    alt_add_custom(mat, k, 0, 0, 1/k, 0, 0);
}

/* Add a scale transformation to matrix. */
void
alt_add_scale(alt_matrix_t *mat, double x, double y)
{
    alt_add_custom(mat, x, 0, 0, y, 0, 0);
}

/* Add a horizontal shear transformation to matrix. */
void
alt_add_hshear(alt_matrix_t *mat, double h)
{
    alt_add_custom(mat, 1, 0, h, 1, 0, 0);
}

/* Add a vertical shear transformation to matrix. */
void
alt_add_vshear(alt_matrix_t *mat, double v)
{
    alt_add_custom(mat, 1, v, 0, 1, 0, 0);
}

/* Add a rotation to matrix. `a` is the angle in radians. */
void
alt_add_rotate(alt_matrix_t *mat, double a)
{
    double c, s;
    c = cos(a); s = sin(a);
    alt_add_custom(mat, c, -s, s, c, 0, 0);
}

/* Add a translation to matrix. */
void
alt_add_translate(alt_matrix_t *mat, double x, double y)
{
    alt_add_custom(mat, 1, 0, 0, 1, x, y);
}

/* Apply affine transformation to all end points. */
void
alt_transform_endpts(alt_array_t *endpts, alt_matrix_t *mat)
{
    alt_endpt_t *endpt;
    double x, y;
    int i;
    for (i = 0; i < (int) endpts->count; i++) {
        endpt = ALT_CAT(alt_endpt_t, endpts, i);
        x = mat->a*endpt->x + mat->c*endpt->y + mat->e;
        y = mat->b*endpt->x + mat->d*endpt->y + mat->f;
        endpt->x = x; endpt->y = y;
    }
}

/* Apply affine transformation to all control points. */
void
alt_transform_ctrpts(alt_array_t *ctrpts, alt_matrix_t *mat)
{
    alt_ctrpt_t *ctrpt;
    double x, y;
    int i;
    for (i = 0; i < (int) ctrpts->count; i++) {
        ctrpt = ALT_CAT(alt_ctrpt_t, ctrpts, i);
        x = mat->a*ctrpt->x + mat->c*ctrpt->y + mat->e;
        y = mat->b*ctrpt->x + mat->d*ctrpt->y + mat->f;
        ctrpt->x = x; ctrpt->y = y;
    }
}
