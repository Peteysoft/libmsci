#include <stdio.h>

#include "time_class.h"
#include "string_petey.h"
#include "peteys_tmpl_lib.h"
#include "simple_temp.h"
#include "dependent_swap.h"
#include "composite_dataset.h"

namespace libpetey {

template<class dt>
long * heapsort_ptr(dt **data, long n);

namespace datasets {

composite_dataset::composite_dataset() {
  n_var=0;
  table_size=1;

  symbol=new string * [1];
  variable=new dataset * [1];
  ndep=new long [1];
  dependencies=new long * [1];

  index=NULL;

  variable[0]=NULL;
  ndep[0]=0;
  symbol[0]=NULL;
  dependencies[0]=NULL;
  nodeleteflag=DELETE_VARS;
}

composite_dataset::composite_dataset(short nodelete) {
  n_var=0;
  table_size=1;

  symbol=new string * [1];
  variable=new dataset * [1];
  ndep=new long [1];
  dependencies=new long * [1];

  index=NULL;

  variable[0]=NULL;
  ndep[0]=0;
  symbol[0]=NULL;
  dependencies[0]=NULL;
  nodeleteflag=nodelete;
}

composite_dataset::~composite_dataset() {
  for (long i=n_var-1; i>=0; i--) {
    if (nodeleteflag == DELETE_VARS && variable[i] != NULL) delete variable[i];
    if (ndep[i] != 0) delete [] dependencies[i];
    if (symbol[i] != NULL) delete symbol[i];
  }
  delete [] symbol;
  delete [] variable;
  if (index != NULL) delete [] index;
  delete [] ndep;
  delete [] dependencies;
}

void composite_dataset::update() {
  if (update_flag==1) return;

  if (index != NULL) delete [] index;
  index=heapsort_ptr(symbol, n_var);

  update_flag=1;
}

long composite_dataset::read(FILE *fileptr) {
  long nread;			//number of bytes read
  long i, j;			//counter
  dstype_type type;		//type of the variable
  long nd;			//number of dependencies
  long *dep_ind;		//the location of the dependencies
  simple_dataset **dep;		//the dependents

  //clean all the old shit up:
  for (long i=n_var-1; i>=0; i--) {
    ///if (variable[i] != NULL) delete variable[i];
    if (ndep[i] != 0) delete [] dependencies[i];
    if (symbol[i] != NULL) delete symbol[i];
  }
  delete [] symbol;
  delete [] variable;
  delete [] ndep;
  delete [] dependencies;

  //read in the number of variables:
  nread=fread(&n_var, 1, sizeof(n_var), fileptr);

  fprintf(datasets_log, "Found %ld variables in file:\n", (int64_t) n_var);

  table_size=n_var;
  symbol=new string *[table_size];
  variable=new dataset *[table_size];

  ndep=new long[table_size];
  dependencies=new long * [table_size];

  //read in the variable names:
  for (i=0;i<n_var;i++) {
    symbol[i]=new string();
    nread+=symbol[i]->read(fileptr);
  }

  update_flag=0;

  update();

  //read in each of the variables in turn:
  for (i=0;i<n_var;i++) {
    //symbol[i]->print();
    //number of dependencies??
    nread+=fread(&nd, 1, sizeof(nd), fileptr);
    fprintf(datasets_log, "%ld dependents found\n", (int64_t) nd);
    ndep[i]=nd;
    if (nd==0) {
      //if there are no dependencies, move straight to the reading:
      nread+=fread(&type, 1, sizeof(type), fileptr);	//read the type
      switch (type) {
        case UNDEF:		//variable is undefined
	  variable[i]=NULL;
	  break;
        case SIMPLE_FLOAT:
          variable[i]=new simple<float>;
          break;
        case SIMPLE_INT32:
          variable[i]=new simple<int32_t>;
          break;
        case SIMPLE_INT64:
          variable[i]=new simple<int64_t>;
          break;
        case SIMPLE_DOUBLE:
          variable[i]=new simple<double>;
          break;
        case SIMPLE_TIME:
          variable[i]=new simple<time_class>;
          break;
        case SIMPLE_STRING:
          variable[i]=new simple<string>;
          break;
        case COMPOSITE:
          variable[i]=new composite_dataset;
      }
    } else {
      //if there are dependencies, find out how many and which ones:
      dep_ind=new long[nd];
      nread+=fread(dep_ind, 1, sizeof(long)*nd, fileptr);
      dependencies[i]=dep_ind;
      dep=new simple_dataset *[nd];
      for (j=0; j<nd; j++) dep[j]=(simple_dataset *) variable[dep_ind[j]];
      //read in the type:
      nread+=fread(&type, 1, sizeof(type), fileptr);
      switch(type) {
        case DEPENDENT_FLOAT:
          variable[i]=new dependent<float>(dep, ndep[i]);
          break;
        case DEPENDENT_INT32:
          variable[i]=new dependent<int32_t>(dep, ndep[i]);
          break;
        case DEPENDENT_DOUBLE:
          variable[i]=new dependent<double>(dep, ndep[i]);
          break;
        case DEPENDENT_FLOAT_S:
          variable[i]=new dependent_swap<float>(dep, ndep[i]);
          break;
        case DEPENDENT_LONG_S:
          variable[i]=new dependent_swap<int32_t>(dep, ndep[i]);
          break;
        case DEPENDENT_DOUBLE_S:
          variable[i]=new dependent_swap<double>(dep, ndep[i]);
          break;
/*        case DEPENDENT_TIME:
          variable[i]=new dependent<time>(dep, ndep[i]);
          break;
        case DEPENDENT_STRING:
          variable[i]=new dependent<string>(dep, ndep[i]);
          break;
*/
      }
      delete [] dep;
    }
    fprintf(datasets_log, "Reading variable of type %d\n", (int32_t) type);

    if (variable[i]!=NULL) nread+=variable[i]->read(fileptr);

    fprintf(datasets_log, "%ld data elements read in\n", (int64_t) variable[i]->nel());

  }

  //print();

  return nread;

}

long composite_dataset::write(FILE *fileptr) {
  long i, j;
  long nwritten;
  dstype_type type;

  fprintf(datasets_log, "Writing out %ld variables to file:\n", (int64_t) n_var);
  //print();

  //write out the number of variables:
  nwritten=fwrite(&n_var, 1, sizeof(n_var), fileptr);

  //write out the variables names:
  for (i=0; i<n_var; i++) {
    nwritten+=symbol[i]->write(fileptr);
  }

  //write out each of the variables in turn:
  for (i=0; i<n_var; i++) {
    //write out the number of dependents:
    nwritten+=fwrite(&ndep[i], 1, sizeof(long), fileptr);
    if (variable[i] == NULL) {
      //write out the type:
      type=UNDEF;
      nwritten+=fwrite(&type, 1, sizeof(type), fileptr);
    } else {
      //write out the dependencies:
      if (ndep[i] != 0) {
        nwritten+=fwrite(dependencies[i], 1, sizeof(long)*ndep[i], fileptr);
      }
      //write out the type:
      type=variable[i]->type_of();
      nwritten+=fwrite(&type, 1, sizeof(type), fileptr);
      //write out the data:
      nwritten+=variable[i]->write(fileptr);
    }
  }

  return nwritten;

}

long composite_dataset::add_var(string name) {
  long i;		//counter
  long loc;		//location of variable
  long index_loc;	//location in index
  string **new_symbol;
  dataset **new_var;
  string **new_index;
  long *new_ndep;
  long **new_dep;

  if (n_var >= table_size) {
    //table size not big enough, new elements must be added:
    //double the current size:
    table_size=2*table_size;

    //the symbol table:
    new_symbol=new string *[table_size];
    new_var=new dataset * [table_size];

    //the dependencies:
    new_ndep=new long[table_size];
    new_dep=new long * [table_size];

    for (i=0;i<n_var;i++) {
      new_symbol[i]=symbol[i];
      new_var[i]=variable[i];
      new_ndep[i]=ndep[i];
      new_dep[i]=dependencies[i];
    }

    for (i=n_var; i<table_size; i++) {
      new_symbol[i]=NULL;
      new_var[i]=NULL;
      new_ndep[i]=0;
      new_dep[i]=NULL;
    }

    delete [] symbol;
    delete [] variable;
    symbol=new_symbol;
    variable=new_var;

    delete [] ndep;
    delete [] dependencies;

    ndep=new_ndep;
    dependencies=new_dep;

  }

  symbol[n_var]=new string(name);
  variable[n_var]=NULL;

  ndep[n_var]=0;

  loc=n_var;
  n_var++;

  update_flag=0;

  return loc;
}


long composite_dataset::search_var(string name, long &loc) {
  long low, mid, high;
  long cmp;

  update();

  low=0;
  high=n_var-1;

  cmp=symbol[index[low]]->cmp(name);

  if (cmp == 0) {
    loc=low;
    return index[low];
  }

  if (cmp > 0) {
    loc=0;
    return -1;
  }

  cmp=symbol[index[high]]->cmp(name);

  if (cmp == 0) {
    loc=high;
    return index[high];
  }

  if (cmp < 0) {
    loc=n_var+1;
    return -1;
  }


  while (high-low > 1) {
    mid=(low+high)/2;
    cmp=symbol[index[mid]]->cmp(name);
    if (cmp < 0) {
      low=mid;
    } else if (cmp > 0) {
      high=mid;
    } else {
      loc=mid;
      return index[mid];
    }
  }

  return -1;

}

long composite_dataset::cvar(long loc, dataset *new_var) {
  long i, j;
  long type;
  long nd, *dep;

  if (loc < 0 || loc >= n_var) return -1;

  type=new_var->type_of();
//  if (dynamic_cast<dependent_dataset *>(new_var)) {
  if (type >= DEPENDENT && type < COMPOSITE) {
    dependent_dataset *dep_var;
    //search for the dependencies and place them in the index:
    dep_var=(dependent_dataset *) new_var;
    nd=dep_var->get_rank();
    fprintf(datasets_log, "Dep. var. rank: %d\n", (int32_t) nd);
    dep=new long[nd];
    for (i=0; i<nd; i++) {
      dep[i]=-1;
      for (j=0; j<n_var; j++) {
        if (dep_var->dependencies[i] == variable[j]) {
	  dep[i]=j;
	  break;
	}
      }
      //if the dependencies aren't already in the table, we cannot place the
      //variable in the table:
      if (dep[i] == -1) {
        char *var_name;
	var_name=(char *) *symbol[loc];
        fprintf(datasets_log, "Insertion failed for symbol %s\n", var_name);
	delete var_name;
        return -1;
      }
    }
    if (ndep[loc] != 0) delete [] dependencies[loc];
    ndep[loc]=nd;
    dependencies[loc]=dep;
  }

  if (variable[loc] != NULL) delete variable[loc];

  variable[loc]=new_var;

  return loc;
}

dataset *composite_dataset::get_var(long loc) {
  if (loc < 0 || loc >= n_var) return NULL;

  return variable[loc];
}

void composite_dataset::print() {
  char *var_name;
  update();

  for (long i=0; i<n_var; i++) {
    if (symbol[index[i]]!= NULL && variable[index[i]] != NULL) {
      var_name=(char *) *symbol[index[i]];
      printf("%d %s %d  ", (int32_t) index[i], var_name, (int32_t) variable[index[i]]->type_of());
      delete [] var_name;
      for (long j=0; j<ndep[index[i]]; j++) printf("%d ", (int32_t) dependencies[index[i]][j]);
      printf("\n");
    }
  }
}

} //end namespace datasets
} //end namespace libpetey

