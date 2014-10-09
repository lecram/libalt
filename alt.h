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

void alt_bbmargin(alt_bbox_t *bb, double mg);
void alt_bound(alt_endpt_t *points, int count, int rounded, alt_bbox_t *bb);

#endif /* _ALT_H */
