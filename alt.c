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
alt_bound(alt_endpt_t *points, int count, int rounded, alt_bbox_t *bb)
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
