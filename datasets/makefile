#
# Copywrite 2004 Peter Mills.  All rights reserved.
#

#########################################################################
########### The user may choose to modify the following macros: #########
#########################################################################

OPT_VER = -g

BASE_PATH=/usr/local

LIB_DIR = $(BASE_PATH)/lib/
INCLUDE_DIR = $(BASE_PATH)/include/

VPATH = $(LIB_DIR):$(INCLUDE_DIR)

CC = g++
CFLAGS = -I$(INCLUDE_DIR) -Wno-deprecated $(OPT_VER)

#########################################################################
################## End modifiable macros. ###############################
#########################################################################

#all: dataset_lib.a dataset_lib.so

all:		libdataset$(OPT_VER).a \
		dataset.h \
		simple_dataset.h \
		simple_temp.h \
		dependent_dataset.h \
		dependent_temp.h \
		dependent_swap.h \
		composite_dataset.h

install:
	cp libdataset$(OPT_VER).a $(LIB_DIR)
#	cp libdataset.so $(LIB_DIR)
	cp composite_dataset.h $(INCLUDE_DIR)
	cp dependent_intc.h $(INCLUDE_DIR)
	cp dataset.h $(INCLUDE_DIR)
	cp simple_dataset.h $(INCLUDE_DIR)
	cp simple_temp.h $(INCLUDE_DIR)
	cp dependent_dataset.h $(INCLUDE_DIR)
	cp dependent_temp.h $(INCLUDE_DIR)
	cp dependent_swap.h $(INCLUDE_DIR)
#	date > date_made_dataset
	
# make the whole library:
libdataset$(OPT_VER).a: composite_dataset$(OPT_VER).o \
		dependent_swap$(OPT_VER).o \
		dependent_temp$(OPT_VER).o \
		dependent_dataset$(OPT_VER).o \
		simple_temp$(OPT_VER).o \
		simple_dataset$(OPT_VER).o \
		heapsort_ptr$(OPT_VER).o \
		dependent_intc$(OPT_VER).o \
		dataset$(OPT_VER).o
	ar -rs 	libdataset$(OPT_VER).a \
		composite_dataset$(OPT_VER).o \
		dependent_intc$(OPT_VER).o \
		heapsort_ptr$(OPT_VER).o \
		dataset$(OPT_VER).o \
		simple_temp$(OPT_VER).o \
		simple_dataset$(OPT_VER).o \
		dependent_dataset$(OPT_VER).o \
		dependent_temp$(OPT_VER).o \
		dependent_swap$(OPT_VER).o

# remove all object files:
clean:
	rm -f *.o
	rm -f libdataset$(OPT_VER).a
	rm -f libdataset$(OPT_VER).so
#	make -C $(CLIB_PATH) clean

# make a dynamic version of the library (not tested):
#dataset_lib.so:	composite_dataset.o \
		dependent_swap.o \
		dependent_temp.o \
		dependent_dataset.o \
		simple_temp.o \
		simple_dataset.o \
		heapsort_ptr.o \
		peteys_tmpl_lib.a
#	$(CC) $(CFLAGS) -shared -o dataset_lib.so \
	 	$(LIB_DIR)libpetey.so \
		heapsort_ptr.o \
		simple_temp.o \
		simple_dataset.o \
		dependent_dataset.o \
		dependent_temp.o \
		dependent_swap.o \
		composite_dataset.o

# template library--contained in another directory:
#libpetey.a:
#	make -C $(CLIB_PATH)

#libpetey.so:
#	make -C $(CLIB_PATH)

# compile the composite dataset class:
composite_dataset$(OPT_VER).o: simple_temp.h dependent_temp.h dataset.h \
			composite_dataset.h composite_dataset.cpp \
			string_petey.h peteys_tmpl_lib.h
	$(CC) $(CFLAGS) -c composite_dataset.cpp -o composite_dataset$(OPT_VER).o

heapsort_ptr$(OPT_VER).o: heapsort_ptr.cc
	$(CC) $(CFLAGS) -c heapsort_ptr.cc -o $@

# compile the dependent dataset template classes:
dependent_swap$(OPT_VER).o: dependent_swap.cpp dependent_swap.h \
			dependent_temp.h dataset.h \
			time_class.h string_petey.h
	$(CC) $(CFLAGS) -c dependent_swap.cpp -o dependent_swap$(OPT_VER).o

# compile the dependent dataset template classes:
dependent_temp$(OPT_VER).o: dependent_temp.cpp dependent_temp.h \
			dependent_dataset.h dataset.h \
			time_class.h string_petey.h
	$(CC) $(CFLAGS) -c dependent_temp.cpp -o dependent_temp$(OPT_VER).o

# don't need as everything has been moved into the header file:
# compile the virtual dependent dataset base class:
dependent_dataset$(OPT_VER).o: dependent_dataset.cpp dependent_dataset.h dataset.h \
			time_class.h string_petey.h
	$(CC) $(CFLAGS) -c dependent_dataset.cpp -o dependent_dataset$(OPT_VER).o

# compile the simple dataset template classes:
simple_temp$(OPT_VER).o: simple_temp.cpp simple_temp.h simple_dataset.h \
			dataset.h dependent_dataset.h \
			time_class.h string_petey.h peteys_tmpl_lib.h
	$(CC) $(CFLAGS) -c simple_temp.cpp -o simple_temp$(OPT_VER).o

# compile the virtual simple base class:
simple_dataset$(OPT_VER).o: simple_dataset.cpp simple_dataset.h dataset.h \
			dependent_dataset.h
	$(CC) $(CFLAGS) -c simple_dataset.cpp -o simple_dataset$(OPT_VER).o

dependent_intc$(OPT_VER).o: dependent_intc.cc dependent_intc.h \
			dataset.h simple_dataset.h simple_temp.h \
			dependent_dataset.h
	$(CC) $(CFLAGS) -c dependent_intc.cc -o $@

# all this has now been inlined and so is contained in the header file:
# all the inlined functions have now been "un-inlined..."
dataset$(OPT_VER).o: dataset.cpp dataset.h
	$(CC) $(CFLAGS) -c dataset.cpp -o dataset$(OPT_VER).o

#peteys_tmpl_lib.h string_petey.h time_class.h:
#	Make -C $(CLIB_PATH)

