#
# This software is released under the following terms:
#
# 1. No commercial use.
# 2. Copies and derivative works are free to use and modify.
# 3. Attribution must be given to all contributors of both original and derivative works.
#
# Authors:
#
# 2017-07-16 Peter Mills: added license information 
#
# libpetey level 1 Makefile
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

#executable extension:
EXE_EXT = .exe

# paths -- library, include, binary executables
BASE_PATH = /home/lenovo
#BASE_PATH = /home/Peter
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
GSL_LIB = /usr/local/lib

# C++ compiler options:
#CFLAGS = $(OPT) -Wno-deprecated -I$(INCLUDE_PATH) -I$(GSL_INCLUDE) -g
CFLAGS = $(OPT) -I$(INCLUDE_PATH) -I$(GSL_INCLUDE)

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
		symbol_table.h \
		full_util.h \
		kselect.h \
		erfinv.h \
		roots_mins.h \
		solve_lode.h \
		constrained.h \
		gsl_util.h \
		rk_dumb_ts.h \
		peteys_no_templates.h \
		interpol_coeff.h

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
		grammar_t$(OPT).o \
		erfinv$(OPT).o \
		roots_mins$(OPT).o \
		solve_lode$(OPT).o \
		constrained$(OPT).o \
		gsl_util$(OPT).o \
		rk_dumb_ts$(OPT).o \
		peteys_no_templates$(OPT).o \
		interpol_coeff$(OPT).o

#list of executables:
EXE = ignore_code
EXE2=$(addsuffix $(EXE_EXT), $(EXE))

#list of executables in the test suite:
TEST = sorter \
		random \
		test_parse_opts \
		lp_test_suite \
		test_libpetey
		
TEST2=$(addsuffix $(EXE_EXT), $(TEST))

all: libpetey$(OPT).a $(EXE) date_calc 

petey: libpetey$(OPT).a

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

clean_all:
	make OPT=-g clean
	make OPT=-pg clean
	make OPT=-O2 clean

install: libpetey$(OPT).a $(EXE) date_calc
	install libpetey$(OPT).a $(LIB_PATH)
	install $(INCLUDE) $(INCLUDE_PATH)
	install $(EXE) $(BIN_PATH)
	install date_calc $(BIN_PATH)
	install date_calc.1 $(MANPATH)/man1

install_allopt: libpetey-g.a libpetey-O2.a libpetey-pg.a
	make install OPT=-g
	make install OPT=-pg
	make install OPT=-O2

test: $(TEST2)
	./test_libpetey
	./lp_test_suite

# build the static library:
libpetey$(OPT).a:	$(OBJECT)
	ar -rs libpetey$(OPT).a $(OBJECT)

# this part is built rather oddly (*** time for a refactor? ***)
peteys_tmpl_lib.cc: heapsort_tmpl.cc bin_search.cc rk_dumb_tmpl.cc treesort.cpp
	touch peteys_tmpl_lib.cc

$(OBJECT): %$(OPT).o: %.cc $(INCLUDE)
	$(CPP) $(CFLAGS) -c $< -o $@


$(EXE2) $(TEST2): %$(EXE_EXT): %.cc libpetey$(OPT).a
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


