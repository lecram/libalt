#include "alt.h"

void
alt_bbmargin(alt_bbox_t *bb, double mg)
{
    bb->x0 -= mg;
    bb->y0 -= mg;
    bb->x1 += mg;
    bb->y1 += mg;
}
