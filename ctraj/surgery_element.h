#include "tcoord_defs.h"

extern int comp_se(surgery_element *se1, surgery_element *se2);

class surgery_element {
  friend class surgery_obj;
  friend comp_se;

  protected:
    float x;
    float y;
    short hemi;
    char touch_flag;

    surgery_element *next_adj;
    surgery_element *next_sort;
    surgery_element *prev_sort;

    static float xbin;
    static float ybin;
    static float x0;
    static float y0;

  public:
    surgery_element(float xx, float yy);
    ~surgery_element;
    void set_binsize(float xb, float yb, float xz, float yz);
    char touch_on();
    void touch_off();
    
    long sort();		//moves itself up or down in the list
				//to maintain sorted order
				//returns number of places moved
    long fix(float mind);	//adds new points between itself
				//and adjacent in proportion to
				//d/dmin where d is the distance
				//returns the number of new points

};

inline char surgery_element::touch_on() {
  char old;
  old=touch_flag;
  touch_flag=1;
  return old;
}

inline char surgery_element::touch_off() {
  touch_flag=0;
}


