#
# libagf level 2 Makefile
#

#---------------------------------------------------------------------#
#            The user should modify the following macros              #
#            to suit his or her system:                               #
#---------------------------------------------------------------------#

OPT_VER=-g

#EXE_EXT=.exe

CC = g++

BASE_PATH = /home/peter

LIB_PATH = $(BASE_PATH)/lib
BIN_PATH = $(BASE_PATH)/bin
INCLUDE_PATH = $(BASE_PATH)/include

#GSL_LIB = $(LIB_PATH)
#GSL_INCLUDE = $(INCLUDE_PATH)

GSL_LIB = /usr/local/lib
GSL_INCLUDE = /usr/local/include

FORTRAN_RUNTIME = gfortran

LD_OPTIONS = -L$(GSL_LIB) -L$(LIB_PATH) -L. \
		-lpetey$(OPT_VER) -lgsl -lgslcblas

CFLAGS = $(OPT_VER) -I$(GSL_INCLUDE) -I$(INCLUDE_PATH) -g

#---------------------------------------------------------------------#
#                      End modifiable macros                          #
#                                                                     #
#---------------------------------------------------------------------#

LIBNAME = libballoon$(OPT_VER).a

LIBPETEY = libpetey$(OPT_VER).a

VPATH = $(GSL_LIB):$(GSL_INCLUDE):$(LIB_PATH)

EXECUTABLES =	agf


OBJECT_FILES = 	balloon_param$(OPT_VER).o \
		knn_param$(OPT_VER).o \
		agf_balloon$(OPT_VER).o 

INCLUDE_FILES = balloon_param.h \
		knn_param.h \
		agf_balloon.h 

EXECUTABLES2=$(addsuffix $(OPT_VER)$(EXE_EXT), $(EXECUTABLES))
EXECUTABLES3=$(addprefix $(BIN_PATH)/, $(addsuffix $(EXE_EXT), $(EXECUTABLES)))

INCLUDE_FILES2=$(addprefix $(INCLUDE_PATH)/, $(INCLUDE_FILES))

all: $(LIBNAME) $(EXECUTABLES2)

install: $(EXECUTABLES3)
	make install_include
	make install_lib
	install validate_pdf.sh $(BIN_PATH)
	install c_borders.sh $(BIN_PATH)

install_lib: $(LIBNAME)
	install $(LIBNAME) $(LIB_PATH)

install_include: $(INCLUDE_PATH)/$(LIBRARY_INCLUDE)

$(INCLUDE_PATH)/$(LIBRARY_INCLUDE): $(LIBRARY_INCLUDE)
	install $(INCLUDE_FILES) $(INCLUDE_PATH)
	install $(LIBRARY_INCLUDE) $(INCLUDE_PATH)

test: $(BIN_PATH)/agf_test_kernel$(EXE_EXT)
	agf_testsuite

$(EXECUTABLES3) $(OLD3): $(BIN_PATH)/%$(EXE_EXT): %$(OPT_VER)$(EXE_EXT)
	install $< $@

clean:
	rm -f $(OBJECT_FILES) 
	rm -f $(EXECUTABLES2) 
	rm -f $(LIBNAME)

$(LIBNAME): $(OBJECT_FILES) Makefile
	ar -rs $(LIBNAME) $(OBJECT_FILES)

$(OBJECT_FILES): %$(OPT_VER).o: %.cc $(LIBRARY_INCLUDE)
	$(CC) $(CFLAGS) -c $< -o $@


