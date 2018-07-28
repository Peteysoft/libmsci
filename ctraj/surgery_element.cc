#include "tcoord_defs.h"
#include "surgery_element.h"

inline int comp_se(surgery_element *se1, surgery_element *se2) {
  long xint1, xint2, yint1, yint2;
  float x2, y2;
  x2=se2->x;
  y2=se2->y;
  if se1->hemi != se2->hemi) tcoord_N2S(x2, y2);
  xint1=(long) ((se1->x-se1->x0)/se1->xbin);
  yint1=(long) ((se1->y-se1->y0)/se1->ybin);
  xint2=(long) ((x2-se2->x0)/se2->xbin);
  yint2=(long) ((y2-se2->y0)/se2->ybin);
  if (xint1>xint2) return -1;
  if (xint1<xint2) return 1;
  if (yint1>yint2) return -1;
  if (yint1<yint2) return 1;
  return 0;
}

surgery_element::~surgery_element() {
  prev_sort->next_sort=next_sort;
  next_sort->prev_sort=prev_sort;
}

long surgery_element::sort() {
  long movement;
  int comp;
  surgery_element *prev;
  surgery_element *next;

  movement=0;
  //if this element is lexically less than the previous, swap the two
  //keep doing this until it is not the case...
  comp=comp_se(this, prev_sort);
  while (comp < 0) {
    next=next_sort;
    prev=prev_sort;
    prev_sort=prev->prev_sort;
    next_sort=prev;
    prev->next_sort=next;
    prev->prev_sort->next_sort=this;
    prev->prev_sort=this;
    next->prev_sort=prev;
    movement--;
    comp=comp_se(this, prev_sort);
  }
  while (comp > 0) {
    next=next_sort;
    prev=prev_sort;
    prev_sort=next;
    next_sort=next->next_sort;
    next->next_sort->prev_sort=this;
    prev->next_sort=next;
    next->prev_sort=prev;
    next->next_sort=this;
    movement++;
    comp=comp_se(this, prev_sort);
  }

  return movement;
}

long surgery_element::fix(float dmin) {
  float x2, y2;
  float dx, dy;
  float s;
  long nnew;
  float sint;
  surgery_element *new_se;

  x2=next_adj->x;
  y2=next_adj->y;
  if (hemi != next_adj->hemi) tcoord_N2S(x2, y2);
  s=tcoord_ds2(x, y, x2, y2);
  dx=x2-x;
  dy=y2-y;

  if (s>dmin) {
    nnew=(long) (s/dmin);
    for (long i=0; i<nnew; i++) {
      //insert the interpolated point into the contour:
      new_se=new surgery_element(x+dx*(float) i/(float) nnew, y+dy*(float) i/(float) nnew);
      new_se->next_adj=next_adj;
      new_se->prev_adj=this;
      next_adj->prev_adj=new_se;
      next_adj=new_se;

      //sort it in the lexical ordering:
      if (i<nnew/2) {
        if (comp_se(this, new_se) < 1) {
          new_se->next_sort=this;
          new_se->prev_sort=prev_sort;
          prev_sort->next_sort=new_se;
          prev_sort=new_se;
        } else {
          new_se->next_sort=next_sort;
          next_sort->prev_sort=new_se;
          new_se->prev_sort=this;
          this->next_sort=new_se;
        }
      } else {
        if (comp_se(next_adj, new_se) < 1) {
          new_se->next_sort=next_adj;
          new_se->prev_sort=next_adj->prev_sort;
          next_adj->prev_sort->next_sort=new_se;
          next_adj->prev_sort=new_se;
        } else {
          new_se->next_sort=next_adj->next_sort;
          next_adj->next_sort->prev_sort=new_se;
          new_se->prev_sort=next_adj;
          next_adj->next_sort=new_se;
        }
      }
      new_se->sort();
    }
  }

  return nnew;

}


