#########################################################################
########### The user may need to modify the following macros: ###########
#########################################################################

# optimization version:
OPT_VER = -g
EXE_EXT=.exe

# compilers:
CC = g++
F77 = gfortran

# paths for binaries and include files:
BASE_PATH = /home/lenovo
#BASE_PATH=/mnt/sda1/home2/pete
LIB_PATH = $(BASE_PATH)/lib
INCLUDE_PATH = $(BASE_PATH)/include
BIN_PATH = $(BASE_PATH)/bin
MANPATH = $(BASE_PATH)/man

# ARPACK directory:
ARPATH = /usr/local/lib

# library names:
FORTRAN_RUNTIME=gfortran
LIBARPACK=arpack_x86
LIBLEX=fl

# location of gsl library (should really have to link it but we do anyway...):
GSL_LIB = /use/local

# compiler directives:
CFLAGS = -Wno-deprecated -I$(INCLUDE_PATH) $(OPT_VER) -D FLEADING_UNDERSCORES= \
		-D FTRAILING_UNDERSCORES=_
FFLAGS = -I$(INCLUDE_PATH) $(OPT_VER) #-fno-underscoring
LD_FLAGS = $(LIBNAME) -L$(LIB_PATH) -lpetey$(OPT_VER) -L$(ARPATH) -l$(LIBARPACK) -l$(FORTRAN_RUNTIME) -L$(GSL_LIB) -lgsl -lgslcblas

# lex:
LEX=flex -I
# yacc:
YY=bison

# location of readline libary:
READLINE_LIB=$(LIB_PATH)

# linker options for sparse_calc:
INT_UTIL_LDF = -L$(READLINE_LIB) -l$(LIBLEX) -lreadline -lncurses


#########################################################################
################## End modifiable macros. ###############################
#########################################################################

VPATH = $(LIB_PATH)

LIBNAME = libsparse$(OPT_VER).a

INCLUDE_FILES = matrix_base.h \
		full_matrix.h \
		sparse_element.h \
		sparse.h \
		sparse_array.h \
		av.h \
		debug.h \
		sparse_solve.h 

OBJECT_FILES = av$(OPT_VER).o \
	       conj_grad$(OPT_VER).o \
	       sparse$(OPT_VER).o \
	       sparse_cnl$(OPT_VER).o \
	       full_matrix$(OPT_VER).o \
	       sparse_array$(OPT_VER).o \
	       matrix_base$(OPT_VER).o \
	       sparse_invert$(OPT_VER).o \
	       sparse_arpack_solver$(OPT_VER).o

FORTRAN_OBJ = arevd$(OPT_VER).o \
	      arsvd$(OPT_VER).o

EXECUTABLES = ev_sparse_array \
	      sparse_mat_prod \
	      sparse_to_full_prod \
	      sparse_vect_prod \
	      svd_sparse_array

SPARSE_CALCULABLES = sparse_calc \
		sparse_calc.tab$(OPT_VER).o \
		lex.yy$(OPT_VER).o \
		sparse_calc_defs$(OPT_VER).o

EXECUTABLES2=$(addsuffix $(OPT_VER)$(EXE_EXT), $(EXECUTABLES))
EXECUTABLES3=$(addprefix $(BIN_PATH)/, $(addsuffix $(EXE_EXT), $(EXECUTABLES)))

%INCLUDE_FILES2=$(addprefix $(INCLUDE_PATH)/, $(INCLUDE_FILES))

all: $(LIBNAME) $(EXECUTABLES2) sparse_calc

install: $(EXECUTABLES3)
	make install_include
	make install_lib
	install sparse_calc$(EXE_EXT) $(BIN_PATH)
	install sparse_calc.1 $(MANPATH)/man1

install_lib: $(LIBNAME)
	install $(LIBNAME) $(LIB_PATH)

install_include: 
	install $(INCLUDE_FILES) $(INCLUDE_PATH)

$(EXECUTABLES3): $(BIN_PATH)/%$(EXE_EXT): %$(OPT_VER)$(EXE_EXT)
	install $< $@

clean:
	rm -f $(OBJECT_FILES) 
	rm -f $(EXECUTABLES2) 
	rm -f $(LIBNAME)
	rm -f $(SPARSE_CALCULABLES)
	rm -f $(FORTRAN_OBJ)
	rm -f sparse_calc.tab.*
	rm -f lex.yy.c

$(EXECUTABLES2): %$(OPT_VER)$(EXE_EXT): %.cc $(LIBNAME) $(LIBPETEY)
	$(CC) $(CFLAGS) -o $@ $< $(LD_FLAGS)

$(LIBNAME): $(OBJECT_FILES) $(FORTRAN_OBJ) makefile
	ar -rs $(LIBNAME) $(OBJECT_FILES) $(FORTRAN_OBJ)

$(OBJECT_FILES): %$(OPT_VER).o: %.cc $(INCLUDE_FILES)
	$(CC) $(CFLAGS) -c $< -o $@

$(FORTRAN_OBJ): %$(OPT_VER).o: %.f av.h
	$(F77) -c $(FFLAGS) $< -o $@


clean_all:
	make OPT_VER= clean
	make clean
	make OPT_VER=-g clean
	make OPT_VER=-O2 clean
	make OPT_VER=-pg clean

sparse_calc: sparse_calc.tab$(OPT_VER).o lex.yy$(OPT_VER).o sparse_calc_defs$(OPT_VER).o libsparse$(OPT_VER).a
	$(CC) $(CFLAGS) -o sparse_calc sparse_calc.tab$(OPT_VER).o \
		lex.yy$(OPT_VER).o sparse_calc_defs$(OPT_VER).o \
		$(LD_FLAGS) $(INT_UTIL_LDF)
#		-l$(LIBLEX) -lreadline $(LD_FLAGS) -lncurses

lex.yy$(OPT_VER).o: lex.yy.c sparse_calc.tab.h
	$(CC) -c $(CFLAGS) lex.yy.c -o $@

sparse_calc.tab$(OPT_VER).o: sparse_calc.tab.c
	$(CC) -c $(CFLAGS) sparse_calc.tab.c -o $@

sparse_calc.tab.c sparse_calc.tab.h: sparse_calc.y sparse_calc_defs.h
	$(YY) -v -d sparse_calc.y

lex.yy.c: sparse_lex.l sparse_calc_defs.h
	$(LEX) sparse_lex.l

sparse_calc_defs$(OPT_VER).o: sparse_calc_defs.cc sparse_calc_defs.h \
			av.h sparse_calc_defs.h libpetey$(OPT_VER).a
	$(CC) -c $(CFLAGS) sparse_calc_defs.cc -o $@
