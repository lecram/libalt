#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>

#include "alt.h"

/* Create a new array. Return NULL if there is not enough memory. */
AltArray *
alt_new_array(size_t item_size, unsigned int init_bulk)
{
    AltArray *array;
    array = malloc(sizeof(*array));
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
alt_resize_array(AltArray *array, unsigned int min_bulk)
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
alt_push(AltArray *array, void *item)
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
alt_pop(AltArray *array, void *item)
{
    array->count--;
    if (item != NULL)
        memcpy(item, ALT_AT(array, array->count), array->size);
    if (array->count < array->bulk >> 1)
        alt_resize_array(array, array->count);
}

/* Sort array according to `comp` (see stdlib's qsort()). */
void
alt_sort(AltArray *array, int (*comp)(const void *, const void *))
{
    qsort(array->items, (size_t) array->count, array->size, comp);
}

/* Delete `array` and its content from memory. */
void
alt_del_array(AltArray **array)
{
    free((*array)->items);
    free(*array);
    *array = NULL;
}

AltArray *
alt_box_array(unsigned int item_count, size_t item_size, void *items)
{
    AltArray *array;
    array = malloc(sizeof(*array));
    if (array == NULL) return NULL;
    array->bulk = item_count;
    array->count = item_count;
    array->size = item_size;
    array->items = items;
    return array;
}

void *
alt_unbox_array(AltArray **array)
{
    void *items;
    items = (*array)->items;
    free(*array);
    *array = NULL;
    return items;
}

/* Intersection comparison for qsort()-like interfaces. */
int
alt_comp_cross(const void *a, const void *b)
{
    const AltCross *arg1 = a;
    const AltCross *arg2 = b;
    return arg1->dist - arg2->dist;
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
AltImage *
alt_new_image(int width, int height)
{
    AltImage *image;
    int i;
    image = malloc(sizeof(*image));
    if (image == NULL) return NULL;
    image->width  = width;
    image->height = height;
    image->x0 = image->y0 = INT_MAX;
    image->x1 = image->y1 = INT_MIN;
    image->diameter = 1;
    image->hori = malloc(image->height * sizeof(*image->hori));
    image->vert = malloc(image->width  * sizeof(*image->vert));
    image->extr = malloc(image->height * sizeof(*image->extr));
    if (image->hori == NULL || image->vert == NULL || image->extr == NULL)
        return NULL;
    for (i = 0; i < image->width; i++) {
        image->vert[i] = alt_new_array(sizeof(AltCross), ALT_INIT_BULK);
        if (image->vert[i] == NULL) return NULL;
    }
    for (i = 0; i < image->height; i++) {
        image->hori[i] = alt_new_array(sizeof(AltCross), ALT_INIT_BULK);
        image->extr[i] = alt_new_array(sizeof(AltCross), ALT_INIT_BULK);
        if (image->hori[i] == NULL || image->extr[i] == NULL) return NULL;
    }
    image->data = calloc(4*width*height, sizeof(*image->data));
    if (image->data == NULL) return NULL;
    return image;
}

/* Load image from PAM file. Return NULL on failure. */
AltImage *
alt_open_pam(const char *fname)
{
    char buf[16];
    int depth, maxval;
    AltImage *image;
    int width = 0;
    int height = 0;
    FILE *pam = fopen(fname, "r");
    if (pam == NULL) return NULL;
    fscanf(pam, "%s\n", buf);
    if (strcmp(buf, "P7")) return NULL;
    fscanf(pam, "WIDTH %d\n", &width);
    fscanf(pam, "HEIGHT %d\n", &height);
    fscanf(pam, "DEPTH %d\n", &depth);
    fscanf(pam, "MAXVAL %d\n", &maxval);
    if (width == 0 || height == 0 || depth != 4 || maxval != 255) return NULL;
    fscanf(pam, "TUPLTYPE %s\n", buf);
    if (strcmp(buf, "RGB_ALPHA")) return NULL;
    fscanf(pam, "%s\n", buf);
    if (strcmp(buf, "ENDHDR")) return NULL;
    image = alt_new_image(width, height);
    fread(image->data, sizeof(uint32_t), image->width*image->height, pam);
    fclose(pam);
    return image;
}

/* Set pen diameter (line width). */
double
alt_get_diameter(AltImage *image)
{
    return image->diameter;
}

/* Set pen diameter (line width). */
void
alt_set_diameter(AltImage *image, double diameter)
{
    image->diameter = diameter;
}

/* Get pixel color at (x, y). */
void
alt_get_pixel(AltImage *image, int x, int y,
              uint8_t *r, uint8_t *g, uint8_t *b, uint8_t *a)
{
    int i;
    y = image->height - y - 1;
    i = 4*(y*image->width+x);
    *r = image->data[i];
    *g = image->data[i+1];
    *b = image->data[i+2];
    *a = image->data[i+3];
}

/* Set pixel color at (x, y). */
void
alt_set_pixel(AltImage *image, int x, int y,
              uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
    int i;
    y = image->height - y - 1;
    i = 4*(y*image->width+x);
    image->data[i]   = r;
    image->data[i+1] = g;
    image->data[i+2] = b;
    image->data[i+3] = a;
}

/* Fill the entire image with `color`. */
void
alt_clear(AltImage *image, uint32_t color)
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
alt_blend(AltImage *image, int x, int y,
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
alt_save_pam(AltImage *image, const char *fname)
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
alt_del_image(AltImage **image)
{
    int i;
    for (i = 0; i < (*image)->width; i++) {
        alt_del_array(&(*image)->vert[i]);
    }
    for (i = 0; i < (*image)->height; i++) {
        alt_del_array(&(*image)->hori[i]);
        alt_del_array(&(*image)->extr[i]);
    }
    free((*image)->hori);
    free((*image)->vert);
    free((*image)->extr);
    free((*image)->data);
    free(*image);
    *image = NULL;
}

/* Scan segment pa-pb and add intersections to `image`. */
void
alt_scan(AltImage *image, AltEndPt *pa, AltEndPt *pb)
{
    AltArray *scanline;
    AltCross cross;
    double radius;
    double hxa, hya, hxb, hyb;
    double vxa, vya, vxb, vyb;
    double ex0, ey0, ex1, ey1;
    double slope;
    int hsign, vsign;
    int sx, sy;
    int y;
    radius = image->diameter / 2;
    /*  Define a segment from (vxa, vya) to (vxb, vyb) that coincides with 
     * pa-pb, but is guaranteed to be x-ascending.
     */
    if (pa->x < pb->x) {
        hsign = +1;
        vxa = pa->x; vya = pa->y; vxb = pb->x; vyb = pb->y;
    } else {
        hsign = -1;
        vxa = pb->x; vya = pb->y; vxb = pa->x; vyb = pa->y;
    }
    /*  Define a segment from (hxa, hya) to (hxb, hyb) that coincides with 
     * pa-pb, but is guaranteed to be y-ascending.
     */
    if (pa->y < pb->y) {
        vsign = +1;
        hxa = pa->x; hya = pa->y; hxb = pb->x; hyb = pb->y;
    } else {
        vsign = -1;
        hxa = pb->x; hya = pb->y; hxb = pa->x; hyb = pa->y;
    }
    vxa = round(vxa); vxb = round(vxb);
    hya = round(hya); hyb = round(hyb);
    if (vxa < image->x0) image->x0 = vxa;
    if (vxb > image->x1) image->x1 = vxb;
    if (hya < image->y0) image->y0 = hya;
    if (hyb > image->y1) image->y1 = hyb;
    /*  Add orthogonal intersections to extra scan. This is suboptimal for
     * diagonal segments, specially when fabs(slope) ~ 1.
     */
    ex0 = floor(vxa-radius);
    ey0 = ALT_MAX(0, floor(hya-radius));
    ex1 = ceil(vxb+radius)+1;
    ey1 = ALT_MIN(image->height, ceil(hyb+radius)+1);
    cross.dist = ex0; cross.sign = +1;
    for (y = (int) ey0; y < (int) ey1; y++) {
        scanline = image->extr[y];
        alt_push(scanline, &cross);
    }
    cross.dist = ex1; cross.sign = -1;
    for (y = (int) ey0; y < (int) ey1; y++) {
        scanline = image->extr[y];
        alt_push(scanline, &cross);
    }
    if (pa->y == pb->y) {
        /* Horizontal segment. */
        cross.dist = pa->y; cross.sign = hsign;
        vxa = ALT_MAX(0, vxa);
        vxb = ALT_MIN(image->width, vxb);
        while (vxa < vxb) {
            scanline = image->vert[(int) vxa];
            alt_push(scanline, &cross);
            vxa++;
        }
    } else if (pa->x == pb->x) {
        /* Vertical segment. */
        cross.dist = pa->x; cross.sign = vsign;
        hya = ALT_MAX(0, hya);
        hyb = ALT_MIN(image->height, hyb);
        while (hya < hyb) {
            scanline = image->hori[(int) hya];
            alt_push(scanline, &cross);
            hya++;
        }
    } else {
        /* Diagonal segment. */
        sy = vya < vyb ? +1 : -1;
        slope = (vxb-vxa) / (hyb-hya);
        cross.sign = hsign;
        while (vxa < vxb) {
            if (0 <= (int) vxa && vxa < image->width) {
                scanline = image->vert[(int) vxa];
                cross.dist = vya;
                alt_push(scanline, &cross);
            }
            vya += sy / slope;
            vxa++;
        }
        sx = hxa < hxb ? +1 : -1;
        slope = 1 / slope;
        cross.sign = vsign;
        while (hya < hyb) {
            if (0 <= (int) hya && hya < image->height) {
                scanline = image->hori[(int) hya];
                cross.dist = hxa;
                alt_push(scanline, &cross);
            }
            hxa += sx / slope;
            hya++;
        }
    }
}

/* Scan `count` points at `points` and add intersections to `image`. */
void
alt_scan_array(AltImage *image, AltEndPt *points, int count)
{
    int i;
    AltEndPt pa, pb;
    for (i = 0; i < count-1; i++) {
        pa = points[i];
        pb = points[i+1];
        alt_scan(image, &pa, &pb);
    }
}

/* Sort intersections and reduce them according to the running winding number.
 * Only crossings that go to or come from a zero winding number are kept.
 */
void
alt_windredux(AltImage *image)
{
    AltArray **scans[3], *scanline;
    AltCross cross, *pcross;
    int count[3];
    int winda, windb;
    int i, j, k;
    scans[0] = image->vert;
    scans[1] = image->hori;
    scans[2] = image->extr;
    count[0] = image->width;
    count[1] = count[2] = image->height;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < count[i]; j++) {
            alt_sort(scans[i][j], alt_comp_cross);
            scanline = alt_new_array(sizeof(AltCross), scans[i][j]->count);
            cross.dist = -HUGE_VAL;
            alt_push(scanline, &cross);
            winda = windb = 0;
            for (k = 0; k < (int) scans[i][j]->count; k++) {
                pcross = ALT_CAT(AltCross, scans[i][j], k);
                windb += pcross->sign;
                if (!winda || !windb) {
                    if (pcross->dist == cross.dist) {
                        /* Duplicate crossings cancel each other. */
                        alt_pop(scanline, NULL);
                    } else {
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

/* Return the minimum distance from `x` to its closest point on `scanline`.
 * Note that x here is a 1D coordinate that follows the scanline's direction.
 * This function is only used as a helper to alt_dist().
 */
static double
alt_scanrange(AltArray *scanline, double x)
{
    AltCross *pcross;
    pcross = (AltCross *) scanline->items;
    while ((++pcross)->dist <= x);
    return ALT_MIN(pcross->dist - x, x - (pcross-1)->dist);
}

/* Return the minimum distance between (x, y) and the path scanned into
 * `image`. If that path doesn't intersect the disk of radius `r` centered
 * at (x, y), then return `r`.
 * This function implements the L1L2 algorithm and is a helper for alt_draw().
 */
static double
alt_dist(AltImage *image, double x, double y, double r)
{
    double mind, mag;
    double j, d, d1, d2;
    int i, xi, yi;
    xi = x;
    yi = y;
    i = 0;
    mind = r*r;
    mag = r/sqrt(2.0);
    do {
        d1 = i*i;
        j = HUGE_VAL;
        if (xi-i >= 0 && xi-i < image->width)
            j = ALT_MIN(j, alt_scanrange(image->vert[xi-i], y));
        if (xi+i >= 0 && xi+i < image->width)
            j = ALT_MIN(j, alt_scanrange(image->vert[xi+i], y));
        if (yi-i >= 0 && yi-i < image->height)
            j = ALT_MIN(j, alt_scanrange(image->hori[yi-i], x));
        if (yi+i >= 0 && yi+i < image->height)
            j = ALT_MIN(j, alt_scanrange(image->hori[yi+i], x));
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
 */
void
alt_draw(AltImage *image, uint32_t fill, uint32_t strk)
{
    AltArray *esl, *hsl;
    AltCross *ecross, *hcross;
    bool border, inside;
    double hlwp, d, m, aa;
    uint8_t fr, fg, fb, fa;
    uint8_t sr, sg, sb, sa;
    int x0, y0, x1, y1, x, y;
    alt_windredux(image);
    alt_unpack_color(fill, &fr, &fg, &fb, &fa);
    alt_unpack_color(strk, &sr, &sg, &sb, &sa);
    hlwp = image->diameter/2 + 0.5;
    x0 = ALT_MAX(0, image->x0 - round(image->diameter));
    y0 = ALT_MAX(0, image->y0 - round(image->diameter));
    x1 = ALT_MIN(image->width, image->x1 + round(image->diameter));
    y1 = ALT_MIN(image->height, image->y1 + round(image->diameter));
    for (y = y0; y < y1; y++) {
        esl = image->extr[y];
        hsl = image->hori[y];
        border = inside = false;
        ecross = (AltCross *) esl->items;
        hcross = (AltCross *) hsl->items;
        ecross++; hcross++;
        for (x = 0; x < image->width; x++) {
            if (x >= ecross->dist) {
                border = !border;
                ecross++;
            }
            if (x >= hcross->dist) {
                inside = !inside;
                hcross++;
            }
            if (x0 <= x && x < x1) {
                if (border) {
                    d = alt_dist(image, x, y, hlwp);
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
                } else if (inside) {
                    alt_blend(image, x, y, fr, fg, fb, fa);
                }
            }
        }
    }
    for (y = 0; y < image->height; y++)
        image->hori[y]->count = image->extr[y]->count = 0;
    for (x = 0; x < image->width; x++)
        image->vert[x]->count = 0;
    image->x0 = image->y0 = INT_MAX;
    image->x1 = image->y1 = INT_MIN;
}

/* Convert `curve` to polyline and add resulting points to array.
 * `endpts` is an array of AltEndPt.
 */
void
alt_add_curve(AltArray *endpts, AltCurve *curve)
{
    /*  We split a Bézier curve into two subcurves repeatedly until we reach
     * almost-straight curves that are then output as straight segments.
     *  Since the segments must be added in the correct order, we need to do
     * a depth-first search as we create a tree of subcurves.
     */
    AltArray *stack;
    AltCurve subcurve;
    AltEndPt a, b, c;
    AltEndPt d, e, f;
    double h;
    stack = alt_new_array(sizeof(AltCurve), ALT_INIT_BULK);
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
        } else {
            /* Add point to polyline. */
            alt_push(endpts, &c);
        }
    }
    alt_del_array(&stack);
}

/* Convert a sequence of linked Bézier curves to a polyline.
 *  `ctrpts` is an array of AltCtrPt. The boolean field `.on` indicates
 * whether a point is on-curve or off-curve. On-curve points are end points;
 * off-curve points are Bézier control points. If two off-curve points appear
 * next to each other, an on-curve point is implied halfway between them.
 * Return an array of AltEndPt. */
AltArray *
alt_unfold(AltCtrPt *ctrpts, int count)
{
    AltCtrPt *s, *a, *b, m;
    AltArray *endpts;
    AltCurve curve;
    int i;
    m.on = true;
    endpts = alt_new_array(sizeof(AltEndPt), ALT_INIT_BULK);
    s = &ctrpts[0];
    b = &ctrpts[1];
    alt_push(endpts, s);
    if (b->on) {
        alt_push(endpts, b);
        s = b;
    }
    for (i = 2; i < count; i++) {
        a = b;
        b = &ctrpts[i];
        if (a->on) {
            if (b->on) {
                alt_push(endpts, b);
                s = b;
            }
        } else {
            if (b->on) {
                curve.a.x = s->x; curve.a.y = s->y;
                curve.b.x = a->x; curve.b.y = a->y;
                curve.c.x = b->x; curve.c.y = b->y;
                alt_add_curve(endpts, &curve);
                s = b;
            } else {
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
AltCtrPt *
alt_circle(double x, double y, double r)
{
    double a;
    AltCtrPt *ctrpts;
    a = r * (sqrt(2) - 1);
    ctrpts = malloc(10 * sizeof(*ctrpts));
    ctrpts[0].x = x  ; ctrpts[0].y = y+r; ctrpts[0].on = true;
    ctrpts[1].x = x+a; ctrpts[1].y = y+r; ctrpts[1].on = false;
    ctrpts[2].x = x+r; ctrpts[2].y = y+a; ctrpts[2].on = false;
    ctrpts[3].x = x+r; ctrpts[3].y = y-a; ctrpts[3].on = false;
    ctrpts[4].x = x+a; ctrpts[4].y = y-r; ctrpts[4].on = false;
    ctrpts[5].x = x-a; ctrpts[5].y = y-r; ctrpts[5].on = false;
    ctrpts[6].x = x-r; ctrpts[6].y = y-a; ctrpts[6].on = false;
    ctrpts[7].x = x-r; ctrpts[7].y = y+a; ctrpts[7].on = false;
    ctrpts[8].x = x-a; ctrpts[8].y = y+r; ctrpts[8].on = false;
    ctrpts[9].x = x  ; ctrpts[9].y = y+r; ctrpts[9].on = true;
    return ctrpts;
}

/* Reset matrix to identity map. */
void
alt_reset(AltMatrix *mat)
{
    mat->a = 1; mat->c = 0; mat->e = 0;
    mat->b = 0; mat->d = 1; mat->f = 0;
}

/* Add a custom transformation to matrix. */
void
alt_add_custom(AltMatrix *mat, double a, double b,
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
alt_add_squeeze(AltMatrix *mat, double k)
{
    alt_add_custom(mat, k, 0, 0, 1/k, 0, 0);
}

/* Add a scale transformation to matrix. */
void
alt_add_scale(AltMatrix *mat, double x, double y)
{
    alt_add_custom(mat, x, 0, 0, y, 0, 0);
}

/* Add a horizontal shear transformation to matrix. */
void
alt_add_hshear(AltMatrix *mat, double h)
{
    alt_add_custom(mat, 1, 0, h, 1, 0, 0);
}

/* Add a vertical shear transformation to matrix. */
void
alt_add_vshear(AltMatrix *mat, double v)
{
    alt_add_custom(mat, 1, v, 0, 1, 0, 0);
}

/* Add a rotation to matrix. `a` is the angle in radians. */
void
alt_add_rotate(AltMatrix *mat, double a)
{
    double c, s;
    c = cos(a); s = sin(a);
    alt_add_custom(mat, c, -s, s, c, 0, 0);
}

/* Add a translation to matrix. */
void
alt_add_translate(AltMatrix *mat, double x, double y)
{
    alt_add_custom(mat, 1, 0, 0, 1, x, y);
}

/* Apply affine transformation to all end points. */
void
alt_transform_endpts(AltEndPt *endpts, int count, AltMatrix *mat)
{
    AltEndPt *endpt;
    double x, y;
    int i;
    for (i = 0; i < count; i++) {
        endpt = &endpts[i];
        x = mat->a*endpt->x + mat->c*endpt->y + mat->e;
        y = mat->b*endpt->x + mat->d*endpt->y + mat->f;
        endpt->x = x; endpt->y = y;
    }
}

/* Apply affine transformation to all control points. */
void
alt_transform_ctrpts(AltCtrPt *ctrpts, int count, AltMatrix *mat)
{
    AltCtrPt *ctrpt;
    double x, y;
    int i;
    for (i = 0; i < count; i++) {
        ctrpt = &ctrpts[i];
        x = mat->a*ctrpt->x + mat->c*ctrpt->y + mat->e;
        y = mat->b*ctrpt->x + mat->d*ctrpt->y + mat->f;
        ctrpt->x = x; ctrpt->y = y;
    }
}
