typedef int natural;

natural *rr_work;

struct row_fill {
  natural n;
  natural *fill;
};

row_fill *empty_row() {
  row_fill *result;
  result=new row_fill;
  result->n=0;
  result->fill=NULL;
  return result;
}

row_fill *copy_row_fill(row_fill *old) {
  row_fill *result;
  result=new row_fill;
  result->n=old->n;
  result->fill=new natural[old->n];
  for (natural i=0; i<old->n; i++) result->fill[i]=old->fill[i];
  return result;
}

natural insert_row_element(row_fill *row, natural column) {
  natural *new_fill;
  natural extra=0;

  new_fill=new natural[row->n+1];
  if (row->n == 0) {
    row->fill=new natural[1];
    row->fill[0]=column;
    extra=1;
    row->n++;
  }

  if (row->fill[0] > column) {
    new_fill[0]=column;
    row->n++;
    extra=1;
    for (natural i=1; i<row->n; i++) new_fill[i]=row->fill[i-1];
  } else {
    for (natural i=0; row->fill[i]<=column; i++) {
      new_fill[i]=row->fill[i];
    }
    if (row->fill[i]!=column) {
      row->n++; 
      extra=1;
    }
    for (natural i=extra; i<row->n; i++) new_fill[i]=row->fill[i-extra];
  }

  delete [] row->fill;
  row->fill=new_fill;

  return extra;
}

void delete_row_fill(fow_filling *row) {
  if (row-n > 0) delete [] row_fill->fill;
  delete row;
}

row_fill * fill_intersection(row_fill *row1, row_fill *row2) {
  natural i, j;
  natural k;
  row_fill *nrf;
  nrf=new row_fill;
  nrf->n=0;
  i=0; j=0;
  while (i<row1->n && j<row2->n) {
    if (row1->fill[i] < row2->fill[i]) {
      i++;
    } else if (row1->fill[i] > row2->fill[i]) {
      j++;
    } else {
      rr_work[nrf->n]=i;
      nrf->n++;
    }
  }

  nrf->fill=new natural[nrf->n];
  for (i=0; i<nrf->n; i++) nrf->fill[i]=work[i];

  return nrf;
}

row_fill * fill_union(row_fill *row1, row_fill *row2) {
  natural i, j;
  natural k;
  row_fill *nrf;
  nrf=new row_fill;
  nrf->n=0;
  i=0; j=0;
  while (i<row1->n && j<row2->n) {
    if (row1->fill[i] < row2->fill[i]) {
      rr_work[nrf->n]=i;
      i++;
    } else if (row1->fill[i] > row2->fill[i]) {
      rr_work[nrf->n]=j;
      j++;
    } else {
      rr_work[nrf->n]=i;
    }
    nrf->n++;
  }

  nrf->fill=new natural[nrf->n];
  for (i=0; i<nrf->n; i++) nrf->fill[i]=work[i];

  return nrf;
}

row_fill * reduce_two_rows(row_fill *row1, row_fill *row2, natural delete_column) {
  natural i, j;
  natural k1, k2;
  row_fill *nrf;
  nrf=new row_fill;
  nrf->n=0;
  i=0; j=0;
  while (i<row1->n && j<row2->n) {
    if (row1->fill[i] < row2->fill[j]) {
      rr_work[nrf->n]=i;
      i++;
    } else if (row1->fill[i] > row2->fill[j]) {
      rr_work[nrf->n]=j;
      j++;
    } else {
      if (i!=delete_column) {
        rr_work[nrf->n]=i;
        nrf->n++;
      }
    }
  }

  nrf->fill=new natural[nrf->n];
  for (i=0; i<nrf->n; i++) nrf->fill[i]=rr_work[i];

  return nrf;
}

struct matrix_fill {
  natural m, n;
  row_fills *row;
  boolean *orig_flag;		//is the row original? if so delete on clean up
};

matrix_fill *empty_matrix(natural m1, natural n1) {
  matrix_fill *result;
  result=new matrix_fill;
  result->m=m1;
  result->n=n1;
  result->row=new row_fill[m1];
  orig_flag=new boolean[old->m];
  for (natural i=0; i<m1; i++) {
    result->row[i]=empty_row();
    orig_flag[i]=1;
  }
  return result;
}

natural fill_matrix(matrix_fill *mf, natural i, natural j) {
  return insert_row_element(mf->row[i], j);
}

matrix_fill *copy_matrix_fill(matrix_fill *old) {
  matrix_fill *result;
  result=new matrix_fill;
  result->m=old->m;
  result->n=old->n;
  result->row=new row_fill[old->m];
  orig_flag=new boolean[old->m];
  for (natural i=0; i<old->m; i++) {
    result->row[i]=old->row[i];
    orig_flag[i]=0;
  }
  return result;
}

void delete_matrix_fill(matrix_fill *mf) {
  for (natural i=0; i<mf->m; i++) if (orig_flag[i]) delete_row_fill(mf->row[i]);
  delete [] row;
  delete [] orig_flag;
}

natural n_filled(matrix_fill *mf) {
  natural result=0;
  for (natural i=0; i<mf->m; i++) result+=mf->row[i]->n;
  return result;
}

int matrix_is_singular(matrix_fill *mf) {
  for (natural i=0; i<mf->n; i++) rr_work[i]=0;
  for (natural i=0; i<mf->m; i++) {
    if (mf->row[i]->n==0) return 1;
    for (natural j=0; j<mf->row[i]->n; j++) rr_work[mf->row[i]->fill[j]]=1;
  }
  for (natural i=0; i<mf->n; i++) if (rr_work[i]==0) return 0;
  return 1;
}

matrix_fill * row_reduce_step(matrix_fill *mf, 
		natural i1,		//first row is replaced 
		natural i2, 
		natural j) {
  matrix_fill *result;

  result=copy_matrix_fill(mf);
  result->row[i1]=reduce_two_rows(mf->row[i1], mf->row[i2], j);
  result->orig_flag[i1]=1;
  return result;
}

int row_reduce_step(matrix_fill *mf, natural &i0, natural &j0, natural step[3]) {

  row_fill *intcpt;
  row_fill *unn;

  for (i=i0; i<mf->m; i++) {
    for (j=j0; j<i; j++) {
      intcpt=new row_fill;
      unn=new row_fill;
      intcp->fill=new natural[matrix_fill->row[i]->n+matrix_fill->row[j]];
      if (matrix_fill->row[i]->n > matrix_fill->row[j]->n) {
        unn->fill=new natural[matrix_fill->row[i]->n];
      } else {
        unn->fill=new natural[matrix_fill->row[j]->n];
      }
      intcpt->n=0;
      unn->n=0;
      while (i1<row1->n && j1<row2->n) {
        if (row1->fill[i1] < row2->fill[i1]) {
          intcpt->fill[i1]=row1->fill[i1];
          i1++;
          intcpt->n++;
        } else if (row1->fill[i1] > row2->fill[j1]) {
          intcpt->fill[j1]=row1->fill[j1];
          j1++;
          intcpt->n++;
        } else {
          intcpt->fill[j1]=row1->fill[j1];
          unn->fill[j1]=row1->fill[j1];
          i1++;
          j1++;
          intcpt->n++;
          unn->n++;
        }
      }
      if (unn->n==0) continue;
      if (intcpt->n > row1->fill[i]+1 && intcpt->n > row1->fill[i]+1) continue;

       
