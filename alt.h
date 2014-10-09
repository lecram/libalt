#ifndef _ALT_H
#define _ALT_H

typedef struct {
    double x0, y0, x1, y1;
} alt_bbox_t;

void alt_bbmargin(alt_bbox_t *bb, double mg);

#endif /* _ALT_H */
