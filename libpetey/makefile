#
# Copyright 2004 Peter Mills.  All rights reserved.
#
# The make file for building a set of templated libraries of various
# algorithmic constructs.  Chief amongst these are: binary search
# routines, heapsort routines, a class for working with dates and times,
# a Runge-Kutta integrator and a linked list class.
# 

#####################################################################
########### The user must modify the following macros  ##############
########### to suit his or her system....              ##############
#####################################################################

# optimization level:
OPT=-g

# C++ compiler:
CPP = g++

BASE_PATH = /home/lenovo
#BASE_PATH = /mnt/sda1/home2/pete
#BASE_PATH = /home/pmills

LIB_PATH = $(BASE_PATH)/lib
INCLUDE_PATH = $(BASE_PATH)/include
BIN_PATH = $(BASE_PATH)/bin

# man path:
MANPATH=$(BASE_PATH)/man
#MANPATH=/usr/local/man

# GSL include location:
GSL_INCLUDE = /usr/include
# GSL library locations:
GSL_LIB = /usr/lib

# C++ compiler options:
#CFLAGS = $(OPT) -Wno-deprecated -I$(INCLUDE_PATH) -I$(GSL_INCLUDE) -g
CFLAGS = $(OPT) -I$(INCLUDE_PATH) -I$(GSL_INCLUDE) -g

#linker flags (standard executables):
LDFLAGS = -L$(GSL_LIB) -lgsl -lgslcblas

#the following macros are for executables built with lex/yacc:
# readline library:
READLINE_LIB=$(BASE_PATH)/lib

# lex & yacc:
LEX = flex -I
YY  = bison
LIBLEX=fl

# linker options for interactive utilities sparse_calc and date_calc:
INT_UTIL_LDF = -L$(READLINE_LIB) -l$(LIBLEX) -lreadline -lncurses


#----------------------------------------------------------------#
########### The following macros are for libsparse: ##############
#----------------------------------------------------------------#

# Fortran compiler (for libsparse):
F77 = gfortran
#F77 = g77
FORTRAN_RUNTIME=gfortran
#FORTRAN_RUNTIME=g2c

# ARPACK library (for libsparse):
ARPATH = /usr/local/lib
#ARPATH = $(LIB_PATH)
LIBARPACK=arpack_x86

# extension for executables (for libsparse):
EXE_EXT=.exe

FFLAGS = -I$(INCLUDE_PATH) $(OPT) #-fno-underscoring

# leading and trailing underscores:
FLEADING_UNDERSCORES =
FTRAILING_UNDERSCORES=_

#####################################################################
################   End modifiable macros.     #######################
#####################################################################


VPATH = datasets/ sparse/

all: libpetey$(OPT).a libdataset$(OPT).a libsparse$(OPT).a $(EXE) date_calc 

petey: libpetey$(OPT).a

#list of include files:
INCLUDE = peteys_tmpl_lib.h \
	 	time_class.h \
		string_petey.h \
		linked.h \
		tree_tmp.h \
		bit_array.h \
		tree_lg.h \
		tree_lgi.h \
		kextreme.h \
		quicksort.h \
		supernewton.h \
		randomize.h \
		parse_command_opts.h \
		error_codes.h \
		read_ascii_all.h \
		petey_pointer.h \
		vector_s.h \
		vector_e.h \
		symbol_table.h \
		full_util.h \
		kselect.h \
		erfinv.h \
		roots_mins.h \
		solve_lode.h \
		gsl_util.h \
		rk_dumb_ts.h \
		peteys_no_templates.h

#list of object files/modules:
OBJECT =	symbol_table$(OPT).o \
		peteys_tmpl_lib$(OPT).o \
		tree_tmp$(OPT).o \
		time_class$(OPT).o \
		read_ascii_all$(OPT).o \
		linked$(OPT).o \
		bit_array$(OPT).o \
		supernewton$(OPT).o \
		kextreme$(OPT).o \
		quicksort$(OPT).o \
		tree_lg$(OPT).o \
		tree_lgi$(OPT).o \
		parse_command_opts$(OPT).o \
		randomize$(OPT).o \
		string_petey$(OPT).o \
		full_util$(OPT).o \
		kselect$(OPT).o \
		vector_s$(OPT).o \
		grammar_t$(OPT).o \
		erfinv$(OPT).o \
		roots_mins$(OPT).o \
		solve_lode$(OPT).o \
		gsl_util$(OPT).o \
		rk_dumb_ts$(OPT).o \
		peteys_no_templates$(OPT).o

#list of executables:
EXE = ignore_code

#list of executables in the test suite:
TEST = sorter \
		random \
		test_parse_opts \
		lp_test_suite
		

libdataset$(OPT).a: install_petey
	make -C datasets LIB_DIR=$(LIB_PATH) INCLUDE_DIR=$(INCLUDE_PATH) \
		OPT_VER=$(OPT) CC=$(CPP) CFLAGS="$(CFLAGS)"

libsparse$(OPT).a: install_petey
	make -C sparse LIB_PATH=$(LIB_PATH) INCLUDE_PATH=$(INCLUDE_PATH) \
		BIN_PATH=$(BIN_PATH) MANPATH=$(MANPATH) \
		OPT_VER=$(OPT) EXE_EXT=$(EXE_EXT) OPT_VER=$(OPT) CC=$(CPP) \
		LEX=$(LEX) YY=$(YY) INT_UTIL_LDF="$(INT_UTIL_LDF)" \
		CFLAGS="$(CFLAGS) -DFLEADING_UNDERSCORES=$(FLEADING_UNDERSCORES) -DFTRAILING_UNDERSCORES=$(FTRAILING_UNDERSCORES)" \
		ARPATH=$(ARPATH) LIBARPACK=$(LIBARPACK) \
		F77=$(F77) FFLAGS="$(FFLAGS)" FORTRAN_RUNTIME=$(FORTRAN_RUNTIME) \
		GSL_LIB=$(GSL_LIB)

allopt:
	make OPT=-g all
	make OPT=-pg all
	make OPT=-O2 all

clean:
	rm -f libpetey*.a
	rm -f libpetey*.so
	rm -f $(OBJECT)
	rm -f $(EXE)
	rm -f $(TEST)
	rm -f date_calc.tab.*
	rm -f date_calc
	rm -f lex.yy.c
	make -C sparse OPT_VER=$(OPT) clean
	make -C datasets OPT_VER=$(OPT) clean

clean_all:
	make OPT=-g clean
	make OPT=-pg clean
	make OPT=-O2 clean

install_petey: libpetey$(OPT).a ignore_code date_calc
	install libpetey$(OPT).a $(LIB_PATH)
	install $(INCLUDE) $(INCLUDE_PATH)
	install $(EXE) $(BIN_PATH)
	install date_calc $(BIN_PATH)
	install date_calc.1 $(MANPATH)/man1

install: install_petey
	make install -C datasets LIB_DIR=$(LIB_PATH) INCLUDE_DIR=$(INCLUDE_PATH) \
			OPT_VER=$(OPT)
	make install -C sparse LIB_PATH=$(LIB_PATH) INCLUDE_PATH=$(INCLUDE_PATH) \
			BIN_PATH=$(BIN_PATH) MANPATH=$(MANPATH) \
			OPT_VER=$(OPT) EXE_EXT=$(EXE_EXT)
	cp ignore_code $(BIN_PATH)

install_allopt: libpetey-g.a libpetey-O2.a libpetey-pg.a
	make install_petey OPT=-g
	make install_petey OPT=-pg
	make install_petey OPT=-O2
	make install -C datasets LIB_DIR=$(LIB_PATH) INCLUDE_DIR=$(INCLUDE_PATH) OPT_VER=-g
	make install -C sparse LIB_PATH=$(LIB_PATH) INCLUDE_PATH=$(INCLUDE_PATH) \
			BIN_PATH=$(BIN_PATH) MANPATH=$(MANPATH) \
			EXE_EXT=$(EXE_EXT) OPT_VER=-g
	make install -C datasets LIB_DIR=$(LIB_PATH) INCLUDE_DIR=$(INCLUDE_PATH) OPT_VER=-O2
	make install -C sparse LIB_PATH=$(LIB_PATH) INCLUDE_PATH=$(INCLUDE_PATH) \
			BIN_PATH=$(BIN_PATH) MANPATH=$(MANPATH) \
			EXE_EXT=$(EXE_EXT) OPT_VER=-O2
	make install -C datasets LIB_DIR=$(LIB_PATH) INCLUDE_DIR=$(INCLUDE_PATH) OPT_VER=-pg
	make install -C sparse LIB_PATH=$(LIB_PATH) INCLUDE_PATH=$(INCLUDE_PATH) \
			BIN_PATH=$(BIN_PATH) MANPATH=$(MANPATH) \
			EXE_EXT=$(EXE_EXT) OPT_VER=-pg

test: test_libpetey $(TEST)
	./test_libpetey

# build the static library:
libpetey$(OPT).a:	$(OBJECT)
	ar -rs libpetey$(OPT).a $(OBJECT)

# this part is built rather oddly (*** time for a refactor? ***)
peteys_tmpl_lib.cc: heapsort_tmpl.cc bin_search.cc rk_dumb_tmpl.cc treesort.cpp
	touch peteys_tmpl_lib.cc

$(OBJECT): %$(OPT).o: %.cc $(INCLUDE)
	$(CPP) $(CFLAGS) -c $< -o $@

EXE2=$(addsuffix $(EXE_EXT), $(EXE))

$(EXE): %: %.cc $(PETEY)
	$(CPP) $(CFLAGS) -o $@ $< libpetey$(OPT).a $(LDFLAGS)

parse_grammar: parse_grammar.cc libpetey$(OPT).a
	$(CPP) $(CFLAGS) -o $@ parse_grammar.cc libpetey$(OPT).a $(INT_UTIL_LDF)

date_calc: date_calc.tab.c lex.yy$(OPT).o libpetey$(OPT).a
	$(CPP) $(CFLAGS) -o $@ date_calc.tab.c lex.yy$(OPT).o libpetey$(OPT).a \
		$(INT_UTIL_LDF) $(LDFLAGS)

lex.yy$(OPT).o: lex.yy.c date_calc.tab.h
	$(CPP) $(CFLAGS) -c lex.yy.c -o $@

lex.yy.c: date_parse2.l time_class.h date_calc.tab.h
	$(LEX) date_parse2.l

date_calc.tab.c date_calc.tab.h: date_calc.y time_class.h
	$(YY) -v -d date_calc.y


#unit tests:
TEST2=$(addsuffix $(EXE_EXT), $(TEST))

$(TEST): %: %.cc libpetey$(OPT).a libpetey_unit_tests.h
	$(CPP) $(CFLAGS) -o $@ $< libpetey$(OPT).a $(LDFLAGS)

