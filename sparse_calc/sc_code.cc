

sc_error sc_code_interp(symbol_table *st,
		vector_s<int> *vt,
		vector_s<sc_t *> vs,
		vector_s<sc_t *> *stck,
		vector_s<long> *cd) {

  sc_error err;
  int ninstr=cd->size();

  for (int i=0; i<ninstr; i++) {
    int code=(*cd)[i];
    int stckptr=stck->size();

    if (code<SC_N_RESERVED) {
      switch (code) {
        case (SC_CLOBBER_CODE) {
          long id;
          assert((double) *(*vs)[stckptr-1]==2);
          id=(long) ((double) *(*vs)[stckptr-3]);
          vt[id]=(*vs)[stckptr-2];
        }
        case (

