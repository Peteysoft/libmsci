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

VPATH = datasets/ sparse/ libpetey/ libagf/

INSTALL_PETEY = $(LIB_PATH)/libpetey$(OPT).a

all: libpetey$(OPT).a libdataset$(OPT).a libsparse$(OPT).a libagf$(OPT).a

libpetey$(OPT).a:
	make -C libpetey LIB_DIR=$(LIB_PATH) INCLUDE_DIR=$(INCLUDE_PATH) \
		OPT_VER=$(OPT) EXE_EXT=$(EXE_EXT) \
		CC=$(CPP) CFLAGS="$(CFLAGS)" \
		GSL_INCLUDE=$(GSL_INCLUDE) GSL_LIB=$(GSL_LIB) \
		LEX=$(LEX) YY=$(YY) INT_UTIL_LDF="$(INT_UTIL_LDF)"

libdataset$(OPT).a: $(INSTALL_PETEY)
	make -C datasets LIB_DIR=$(LIB_PATH) INCLUDE_DIR=$(INCLUDE_PATH) \
		OPT_VER=$(OPT) CC=$(CPP) CFLAGS="$(CFLAGS)"

libsparse$(OPT).a: $(INSTALL_PETEY)
	make -C sparse LIB_PATH=$(LIB_PATH) INCLUDE_PATH=$(INCLUDE_PATH) \
		BIN_PATH=$(BIN_PATH) MANPATH=$(MANPATH) \
		EXE_EXT=$(EXE_EXT) OPT_VER=$(OPT) \
		CC=$(CPP) CFLAGS="$(CFLAGS)" \
		LEX=$(LEX) YY=$(YY) INT_UTIL_LDF="$(INT_UTIL_LDF)" \
		FLEADING_UNDERSCORES=$(FLEADING_UNDERSCORES) \
		FTRAILING_UNDERSCORES=$(FTRAILING_UNDERSCORES) \
		ARPATH=$(ARPATH) LIBARPACK=$(LIBARPACK) \
		F77=$(F77) FFLAGS="$(FFLAGS)" FORTRAN_RUNTIME=$(FORTRAN_RUNTIME) \
		GSL_LIB=$(GSL_LIB)

libagf$(OPT).a: $(INSTALL_PETEY)
	make -C libagf LIB_DIR=$(LIB_PATH) INCLUDE_DIR=$(INCLUDE_PATH) \
		OPT_VER=$(OPT) EXE_EXT=$(EXE_EXT) \
		CC=$(CPP) CFLAGS="$(CFLAGS)" \
		GSL_INCLUDE=$(GSL_INCLUDE) GSL_LIB=$(GSL_LIB)

allopt:
	make OPT=-g all
	make OPT=-pg all
	make OPT=-O2 all

clean:
	make -C libpetey OPT_VER=$(OPT) clean
	make -C sparse OPT_VER=$(OPT) clean
	make -C datasets OPT_VER=$(OPT) clean

clean_all:
	make OPT=-g clean
	make OPT=-pg clean
	make OPT=-O2 clean

$(INSTALL_PETEY): libpetey$(OPT).a ignore_code$(EXE_EXT) date_calc$(OPT_VER)$(EXE_EXT)
	make install -C libpetey LIB_PATH=$(LIB_PATH) INCLUDE_PATH=$(INCLUDE_PATH) \
			BIN_PATH=$(BIN_PATH) MANPATH=$(MANPATH) \
			OPT_VER=$(OPT) EXE_EXT=$(EXE_EXT)

install: $(INSTALL_PETEY)
	make install -C datasets LIB_DIR=$(LIB_PATH) INCLUDE_DIR=$(INCLUDE_PATH) \
			OPT_VER=$(OPT)
	make install -C sparse LIB_PATH=$(LIB_PATH) INCLUDE_PATH=$(INCLUDE_PATH) \
			BIN_PATH=$(BIN_PATH) MANPATH=$(MANPATH) \
			OPT_VER=$(OPT) EXE_EXT=$(EXE_EXT)
	make install -C libagf LIB_PATH=$(LIB_PATH) INCLUDE_PATH=$(INCLUDE_PATH) \
			BIN_PATH=$(BIN_PATH) \
			OPT_VER=$(OPT) EXE_EXT=$(EXE_EXT)
	cp ignore_code $(BIN_PATH)

install_allopt: libpetey-g.a libpetey-O2.a libpetey-pg.a
	make install OPT=-g
	make install OPT=-pg
	make install OPT=-O2

