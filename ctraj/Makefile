#
# Main contributer: Peter Mills (peteymills@hotmail.com)
#
# Makefile to build a set of executables for modelling trajectory
# advection and tracer transport.
#
# The new models are:
#
# ctraj_model: simple trajectory model.  Can use pre-computed fields,
# 		global fields and analytic fields.
#
# ctraj_tracer: a semi-Lagrangian tracer advection model.  Computes 
#		a mapping from one field to then next and stores them as a
#		series of sparse matrices.
#
# ctraj_contour: a contour advection model.  Uses the curvature criterion for
#		adding and removing points.
#

#_____________________________________________________________________________#
#                                                                             #
#                     The user should modify the following                    #
#                     macros to suit his or her system                        #
#_____________________________________________________________________________#
#                                                                             #

CC = g++

OPT_VER = -g

# paths for binaries and include files:
#BASE_PATH = /home/lenovo
BASE_PATH = /home/peter
#BASE_PATH = /mnt/sdc1/home2/pete/

LIB_PATH = $(BASE_PATH)/lib

INCLUDE_PATH = $(BASE_PATH)/include
BIN_PATH = $(BASE_PATH)/bin

# location of GNU Scientific Library:
GSL = $(BASE_PATH)
GSL_LIB = $(GSL)/lib
GSL_INC = $(GSL)/include

# location of NetCDF libraries:
NETCDF = /usr/local
NCDF_LIB = $(NETCDF)/lib
NCDF_INC = $(NETCDF)/include

# ARPACK library:
ARPATH=/usr/local/lib
LIBARPACK=arpack_x86

# Fortran run-time libary:
#FORTRAN_RUNTIME=g2c
FORTRAN_RUNTIME=gfortran

# compiler directives:
LD_OPTIONS = -L$(GSL_LIB) -L$(LIB_PATH) \
	-ldataset$(OPT_VER) -lsparse$(OPT_VER) -lpetey$(OPT_VER) \
	-lgsl -lgslcblas \
	-L$(ARPATH) -l$(LIBARPACK) -l$(FORTRAN_RUNTIME) -ldl \
	-L$(NCDF_LIB) -lnetcdf_c++ -lnetcdf -lhdf5_hl -lhdf5 -lz
	#-L$(NCDF_LIB) -lnetcdf_c++ -lnetcdf -lz

#LD_OPTIONS = -L$(LIB_PATH) -ldataset$(OPT_VER) -lpetey$(OPT_VER)
CFLAGS = -g -I$(GSL_INC) -I$(NCDF_INC) -I$(INCLUDE_PATH) \
		-Wno-deprecated $(OPT_VER) -D FTRAILING_UNDERSCORES=_

#_____________________________________________________________________________#
#                                                                             #
#                           END modifiable macros                             #
#_____________________________________________________________________________#
#                                                                             #

# library name:
LIBNAME = libctraj.a

# main include file:
LIBRARY_INCLUDE = ctraj.h

# auxilliary libraries:
LIBPETEY = $(LIB_PATH)/libpetey$(OPT_VER).a
LIBSPARSE = $(LIB_PATH)/libsparse$(OPT_VER).a
LIBDATASET = $(LIB_PATH)/libdataset$(OPT_VER).a

VPATH = $(LIB_PATH):$(INCLUDE_PATH)

EXE = ctraj_model \
	ctraj_contour \
	ctraj_tracer \
	tracer_step2_w_renorm \
	proxy_tracer \
	pc_proxy pc_proxy_predict \
	tracer_interpolate vfield_interpolate \
	contour_dimension contour_length eq_lat

UTIL = zonally_symmetric_tracer lla2aeb extract_field \
       		nc2ds ecmwf2ds vdata2ds extract_ncep_field \
		bev2xy xy2bev random_global c0 \
		select_meas correlate_fields \
		split_file 

# object files:
OBJECTS = 	ctraj_defaults.o \
		ctraj_vfield_standard.o \
		ctraj_vfield_2d.o \
		ctraj_vfield_anal.o \
		ctraj_vfield_base.o \
		ctraj_boundary_list.o \
		ctraj_tfield_standard.o \
		ctraj_tfield_diffusion.o \
		ctraj_tfield_nd.o \
		meas_data.o tracer_anal.o \
		metric_base.o az_eq_t.o \
		coordtran.o \
		ncep_util.o \
		ctraj_optargs.o \
		global_metric.o \
		contour_anal.o \
		tcoord_defs.o

#scripts:
SCRIPTS = get_ecmwf.py ecmwf.py ecmwf2ds.sh ncep2ds.sh

# include files:
INCLUDE = ctraj_boundary_list.h \
	ctraj_vfield_base.h \
	ctraj_vfield_standard.h \
	ctraj_vfield_2d.h \
	ctraj_vfield_anal.h \
	ctraj_tfield_base.h \
	ctraj_tfield_standard.h \
	ctraj_tfield_nd.h \
	ctraj_defaults.h \
	metric_base.h az_eq_t.h \
	coordtran.h \
	ncep_util.h tracer_anal.h \
	meas_data.h \
	global_metric.h \
	tcoord_defs.h

all:	$(EXE) $(UTIL) $(LIBNAME)
	make -C scripts/plotting BIN_PATH=$(BIN_PATH)

install:
	install $(EXE) $(BIN_PATH)
	install $(UTIL) $(BIN_PATH)
	install $(LIBNAME) $(LIB_PATH)
	install $(INCLUDE) $(INCLUDE_PATH)
	install $(LIBRARY_INCLUDE) $(INCLUDE_PATH)
	install $(SCRIPTS) $(BIN_PATH)
	make -C scripts/plotting BIN_PATH=$(BIN_PATH) install

#build the main library--this is linked into all the executables...
libctraj.a: $(OBJECTS)
	ar -rs libctraj.a $(OBJECTS)

clean:
	rm -f $(EXE)
	rm -f $(UTIL)
	rm -f $(OBJECTS)
	rm $(LIBNAME)

#compiles the help screens from all the executables into a single text
#helpfile:
HELPFILES=$(addprefix help_, $(EXE) $(UTIL))

helpfile.txt: $(HELPFILES)
	cat $(HELPFILES) > helpfile.txt

$(HELPFILES): help_%: %
	$< -? > $@

# suffix rule to build all the executables:
$(EXE) $(UTIL) $(OLD): %$(EXE_EXT): %.cc $(LIBNAME) $(LIBDATASET) $(LIBPETEY) $(LIBSPARSE) $(LIBRARY_INCLUDE)
	$(CC) $(CFLAGS) -o $@ $< $(LIBNAME) $(LD_OPTIONS)

$(OBJECTS): %.o: %.cc $(LIBRARY_INCLUDE)
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBRARY_INCLUDE): $(INCLUDE)
	touch $(LIBRARY_INCLUDE)

