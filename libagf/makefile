#
# Copyright 2007 Peter Mills.  All rights reserved.
#

#---------------------------------------------------------------------#
#            The user should modify the following macros              #
#            to suit his or her system:                               #
#---------------------------------------------------------------------#

OPT_VER=-g

#EXE_EXT=.exe

CC = g++

#BASE_PATH = /mnt/sda1/home2/pete
BASE_PATH = /home/lenovo

LIB_PATH = $(BASE_PATH)/lib
BIN_PATH = $(BASE_PATH)/bin
INCLUDE_PATH = $(BASE_PATH)/include

#GSL_LIB = $(LIB_PATH)
#GSL_INCLUDE = $(INCLUDE_PATH)

GSL_LIB = /usr/local/lib
GSL_INCLUDE = /usr/local/include

LD_OPTIONS = -L$(GSL_LIB) -L$(LIB_PATH) -L. \
		-lpetey$(OPT_VER) -lgsl -lgslcblas #-L/usr/lib/ncurses -lcurses

CFLAGS = $(OPT_VER) -I$(GSL_INCLUDE) -I$(INCLUDE_PATH) -Wno-deprecated #-O2

#---------------------------------------------------------------------#
#                      End modifiable macros                          #
#                                                                     #
#---------------------------------------------------------------------#

LIBNAME = libagf$(OPT_VER).a
LIBPETEY = libpetey$(OPT_VER).a

TEST_DIR = examples/sample_classes

all:
	make -C src all \
			OPT_VER=$(OPT_VER) CC=$(CC) EXE_EXT=$(EXE_EXT) \
			LIB_PATH=$(LIB_PATH) BIN_PATH=$(BIN_PATH) \
			INCLUDE_PATH=$(INCLUDE_PATH) \
			GSL_LIB=$(GSL_LIB) GSL_INCLUDE=$(GSL_INCLUDE) \
			CFLAGS="$(CFLAGS)" LD_OPTIONS="$(LD_OPTIONS)" \
			LIBPETEY=$(LIBPETEY) LIBNAME=$(LIBNAME)

install:
	make -C src install EXE_EXT=$(EXE_EXT) \
			LIB_PATH=$(LIB_PATH) BIN_PATH=$(BIN_PATH) INCLUDE_PATH=$(INCLUDE_PATH)
#	make -C $(TEST_DIR) install BIN_PATH=$(BIN_PATH) EXE_EXT=$(EXE_EXT)

test: all install
	make -C $(TEST_DIR) all EXE_EXT=$(EXE_EXT) \
			OPT_VER=$(OPT_VER) CC=$(CC) \
			LIB_PATH=$(LIB_PATH) BIN_PATH=$(BIN_PATH) \
			INCLUDE_PATH=$(INCLUDE_PATH) \
			GSL_LIB=$(GSL_LIB) GSL_INCLUDE=$(GSL_INCLUDE) \
			CFLAGS="$(CFLAGS)" LD_FLAGS="$(LD_OPTIONS)" \
			LIBNAME=$(LIBNAME)
	make -C $(TEST_DIR) EXE_EXT=$(EXE_EXT) test 
	make -C $(TEST_DIR) EXE_EXT=$(EXE_EXT) test_pdf 

clean:
	make -C src EXE_EXT=$(EXE_EXT) clean
	make -C $(TEST_DIR) EXE_EXT=$(EXE_EXT) clean
	make -C $(TEST_DIR) clean_test
	make -C $(TEST_DIR) pdf_clean

clean_test:
	make -C $(TEST_DIR) clean_test
	make -C $(TEST_DIR) pdf_clean

