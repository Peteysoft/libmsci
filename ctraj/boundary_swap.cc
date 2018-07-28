#include <string.h>

#include "tcoord_defs.h"
#include "boundary_swap.h"

#define ELSIZE 10

boundary_swap::boundary_swap(char *nfile, char *sfile, char *bfile, long bs) {
  time_class tref1, tref2;

  bsizen=bs;
  tintn=new traj_int_obj(nfile);
  tints=new traj_int_obj(sfile);
  list=new boundary_element;
  list->next=NULL;
  n=0;

  sbase=new char[sizeof(bfile)+1];
  strcpy(sbase, bfile);

  swap=fopen(bfile, "r+");
  nextb=0;
 
  get_next();

  tref1=tintn->get_time(0);
  tref2=tints->get_time(0);

  if (tref2>tref1) tref1=tref2;
  toffs=tints->get_tind(tref1)-tintn->get_tind(tref1);

}

boundary_swap::~boundary_swap() {
  boundary_element *current;
  boundary_element *previous;

  fclose(swap);

  return;

  current=list;
  while (current != NULL) {
    previous=current;
    current=current->next;
    delete previous;
  }

  delete tintn;
  delete tints;
 
}

void boundary_swap::set_parm(double tcoarse, long nfine, float minspc, float maxspc) {

  tstep=tcoarse;
  nrk=nfine;

  min_spac=minspc;
  max_spac=maxspc;

}


time_class boundary_swap::set(interpol_index ind) {
  tind=ind;
  return tintn->get_time(ind);
}

interpol_index boundary_swap::set(time_class t0) {
  tind=tintn->get_tind(t0);
  return tind;
}

ind_type boundary_swap::nt(interpol_index &offset) {
  ind_type nt1, nt2;
  nt1=tintn->nt();
  nt2=tints->nt();

  if (nt2-toffs < nt1) nt1=nt2-toffs;

  offset=toffs;

  return nt1;
}

interpol_index boundary_swap::get_tind() {
  return tintn->get_tind(tind);
}

time_class boundary_swap::get_t() {
  return tintn->get_time(tind);
}

time_class boundary_swap::advance_current() {
  boundary_element *current;
  current=list;

  for (long i=0; i<n; i++) {
    tcoord_fix(current->x, current->y, current->hemi);
    if (current->hemi==1) {
      tintn->integrate((float *) current, (float *) current, tind, tstep/nrk, nrk);
    } else {
      tints->integrate((float *) current, (float *) current, tind+toffs, tstep/nrk, nrk);
    }
    current=current->next;
  }
  tind+=tstep;

  return tintn->get_time(tind);

}

void boundary_swap::fix_current() {
  boundary_element *current;
  boundary_element *next;
  boundary_element *intd;
  float ds;
  long nnew;
  long nins;
  
  current=list;

  nnew=0;

  for (long i=0; i<n-1; i++) {
    next=current->next;
    if (current->hemi != next->hemi) tcoord_N2S(current->x, current->y);
    ds=sqrt(tcoord_ds2(current->x, current->y, next->x, next->y));
    //printf("fix_current: ds=%f\n", ds);
    if (ds < min_spac) {
      current->next=next->next;
      delete next;
      nnew--;
    } else if (ds > max_spac) {
      nins=(long) (ds/max_spac);
      intd=current;
      for (long j=0; j<nins; j++) {
        intd->next=new boundary_element;
        intd=intd->next;
        intd->x=current->x+(next->x-current->x)*(j+1.)/(nins+1.);
        intd->y=current->y+(next->y-current->y)*(j+1.)/(nins+1.);
        intd->hemi=next->hemi;
      }
      intd->next=next;
      current=next;
      nnew+=nins;
    } else {
      current=next;
    }
  }
  printf("fix_current: n=%8d, nnew=%8d\n", n, nnew);

  n+=nnew;
}

void boundary_swap::store_current() {
  time_class tcur;
  boundary_element *current;
  long nextnextb;
  char fname;
  long nrem, bs2;
  long nnew;
  float x, y;
  short hemi;

  current=list;

  tcur=tintn->get_time(tind);

  if (n > bsize) {
    //check where the next block will go:
    fseek(swap, 0, SEEK_END);
    nextnextb=ftell(swap);
    //printf("store_current: nextnextb=%d\n", nextnextb);

    //if the number of elements doesn't exceed twice the block
    //size, then divide them evenly between the two blocks,
    //otherwise, fill the current block and make the next
    //one the size of the remainder...
    nrem=n-bsize;
    if (nrem < bsizen) {
      nnew=n/2;
      nrem=n-nnew;
      bs2=bsizen;
    } else {
      nnew=n;
      bs2=nrem;
    }
    printf("Starting new block: loc=%12d, size=%8d, elements=%8d\n", 
		nextnextb, bs2, nrem);

    //write time, block size, number of elements and location of next block:
    fseek(swap, bstart, SEEK_SET);
    fwrite(&tcur, sizeof(tcur), 1, swap);
    fwrite(&bsize, sizeof(long), 1, swap);
    fwrite(&nnew, sizeof(long), 1, swap);
    fwrite(&nextnextb, sizeof(nextnextb), 1, swap);

    //write the data up to the current block size:
    for (long i=0; i<nnew; i++) {
      fwrite(&current->x, sizeof(float), 1, swap);
      fwrite(&current->y, sizeof(float), 1, swap);
      fwrite(&current->hemi, sizeof(short), 1, swap);
      current=current->next;
    }
    fwrite(&current->x, sizeof(float), 1, swap);
    fwrite(&current->y, sizeof(float), 1, swap);
    fwrite(&current->hemi, sizeof(short), 1, swap);

    //start the new block:
    fseek(swap, 0, SEEK_END);
    fwrite(&tcur, sizeof(tcur), 1, swap);
    fwrite(&bs2, sizeof(long), 1, swap);
    fwrite(&nrem, sizeof(long), 1, swap);
    fwrite(&nextb, sizeof(long), 1, swap);
    
    for (long i=0; i<=nrem; i++) {
      x=current->x;
      y=current->y;
      hemi=current->hemi;
      fwrite(&x, sizeof(float), 1, swap);
      fwrite(&y, sizeof(float), 1, swap);
      fwrite(&hemi, sizeof(short), 1, swap);
      current=current->next;
    }
    //write dummy data to the end of the block:
    for (long i=nrem; i<bs2; i++) {
      fwrite(&x, sizeof(float), 1, swap);
      fwrite(&y, sizeof(float), 1, swap);
      fwrite(&hemi, sizeof(short), 1, swap);
    }
    n=nnew;
  } else {
    //write time, block size and location of next block:
    fseek(swap, bstart, SEEK_SET);
    fwrite(&tcur, sizeof(tcur), 1, swap);
    fwrite(&bsize, sizeof(long), 1, swap);
    fwrite(&n, sizeof(long), 1, swap);
    fwrite(&nextb, sizeof(nextb), 1, swap);
    
    for (long i=0; i<n; i++) {
      fwrite(&current->x, sizeof(float), 1, swap);
      fwrite(&current->y, sizeof(float), 1, swap);
      fwrite(&current->hemi, sizeof(short), 1, swap);
      current=current->next;
    }
  }
}
      
void boundary_swap::get_next () {  
  boundary_element *current, *next;
  long nnew;
  time_class tcur;
  long i, j;
  char tstring[30];

  current=list;

  fseek(swap, nextb, SEEK_SET);
  bstart=nextb;
  fread(&tcur, sizeof(tcur), 1, swap);
  fread(&bsize, sizeof(bsize), 1, swap);
  fread(&nnew, sizeof(nnew), 1, swap);
  fread(&nextb, sizeof(nextb), 1, swap);

  tind=tintn->get_tind(tcur);

  printf("Reading next block: loc=%12d, size=%8d, elements=%8d\n",
		bstart, bsize, nnew);
  //tcur.write_string(tstring);
  //printf("get_next: tcur=%s\n", tstring);

  for (i=0; i<nnew && current->next != NULL; i++) {
    fread(&current->x, sizeof(float), 1, swap);
    fread(&current->y, sizeof(float), 1, swap);
    fread(&current->hemi, sizeof(short), 1, swap);
    current=current->next;
  }

  for (j=i; j<nnew; j++) {
    fread(&current->x, sizeof(float), 1, swap);
    fread(&current->y, sizeof(float), 1, swap);
    fread(&current->hemi, sizeof(short), 1, swap);
    current->next = new boundary_element;
    current=current->next;
    current->next=NULL;
  }

  //clip any extra:
  while (current->next != NULL) {
    next=current->next;
    delete current;
    current=next;
  }

  n=nnew;
}

time_class boundary_swap::advance() {
  double tindn;
  time_class rt;

  tindn=tind+tstep;

  long nblock;

  nblock=0;
  while (tind < tindn) {
    rt=advance_current();
    fix_current();
    store_current();
    get_next();
    nblock++;
  }

  //printf("advance: blocks=%3d, elements=%d, bsize=%d, tind=%lf, tstep=%lf\n", n, bstart, bsize, tind, tstep);

  return rt;
}

long boundary_swap::store(char *fname) {
  long bsref;
  float lon, lat;
  boundary_element *current;
  FILE *fs;
  long ntot;
  int32_t nvar=2;

  bsref=bstart;

  fs=fopen(fname, "w");

  fwrite(&nvar, sizeof(nvar), 1, fs);

  current=list;
  for (long i=0; i<n; i++) {
    tcoord2_2lonlat(current->x, current->y, 1, current->hemi, lon, lat);
    fwrite(&lon, sizeof(float), 1, fs);
    fwrite(&lat, sizeof(float), 1, fs);
    //printf("%f %f\n", lon, lat);
    current=current->next;
  }
  ntot=n;

  get_next();
  while (bstart != bsref) {
    ntot+=n;
    current=list;
    for (long i=0; i<n; i++) {
      tcoord2_2lonlat(current->x, current->y, 1, current->hemi, lon, lat);
      //printf("%d %f %f\n", current->hemi, current->x, current->y);
      fwrite(&lon, sizeof(float), 1, fs);
      fwrite(&lat, sizeof(float), 1, fs);
      current=current->next;
      //printf("%f %f\n", lon, lat);
    }
    get_next();
  }

  fclose(fs);
  return ntot;
}

void boundary_swap::print_current() {
  boundary_element *current;
  current=list;
  for (long i=0; i<n; i++) {
    printf("%f %f %d\n", current->x, current->y, current->hemi);
    current=current->next;
  }

}


