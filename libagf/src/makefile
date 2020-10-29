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

#BASE_PATH = /home/lenovo
BASE_PATH = /home/peter
#BASE_PATH=/mnt/sdc1/home2/pete

LIB_PATH = $(BASE_PATH)/lib
BIN_PATH = $(BASE_PATH)/bin
INCLUDE_PATH = $(BASE_PATH)/include

#GSL_LIB = $(LIB_PATH)
#GSL_INCLUDE = $(INCLUDE_PATH)

GSL_LIB = /usr/local/lib
GSL_INCLUDE = /usr/local/include

FORTRAN_RUNTIME = gfortran

LD_OPTIONS = -L$(GSL_LIB) -L$(LIB_PATH) -L. \
		-lpetey$(OPT_VER) -lgsl -lgslcblas -l$(FORTRAN_RUNTIME) #-L/usr/lib/ncurses -lcurses

#CFLAGS = $(OPT_VER) -I$(GSL_INCLUDE) -I$(INCLUDE_PATH) -Wno-deprecated #-O2 -g
CFLAGS = $(OPT_VER) -I$(GSL_INCLUDE) -I$(INCLUDE_PATH) -g

#---------------------------------------------------------------------#
#                      End modifiable macros                          #
#                                                                     #
#---------------------------------------------------------------------#

LIBSVM = /usr/local/libsvm-3.18

LIBNAME = libagf$(OPT_VER).a

LIBPETEY = libpetey$(OPT_VER).a

VPATH = $(GSL_LIB):$(GSL_INCLUDE):$(LIB_PATH)

EXECUTABLES =	class_borders \
		classify_b \
		classify_c \
		classify_m \
		classify_s \
		svm_accelerate \
		cls_comp_stats \
		pdf_knn \
                pdf_agf \
		pdf_sim \
		pdf_agf_opt \
		nfold \
		cluster_knn \
		agf_precondition \
		agf_preprocess \
		agf_correlate \
		float_to_class \
		class_to_float \
		roc_curve \
		browse_cluster_tree \
		multi_borders \
		C2R \
		lvq2agf \
		agf2ascii \
		svmout2agf \
		svm2agf \
		mbh2mbm \
		agf \
		knn \
		print_control \
		agf_testsuite \
		sort_discrete_vectors \
		search_discrete_vectors \
		optimize_r0 \
		calibrate_borders \
		subsample_special \
		subsample_sv \
		validate_probabilities \
		compute_kernel
# more experimental stuff:
#		optimize_r0 calibrate_p2 calc_joint_entropy \
#		calibrate_p3 class_triangle 

OLD = 		classify_a \
		classify_knn \
		int_agf \
		int_knn \
		test_classify \
		test_classify_b \
		test_knn \

OBJECT_FILES = 	agf_metric2$(OPT_VER).o \
		agf_engine$(OPT_VER).o \
		agf_calc_pdf$(OPT_VER).o \
		agf_calc_pdf_opt$(OPT_VER).o \
		agf_classify$(OPT_VER).o \
		adgaf$(OPT_VER).o \
		sample_class_borders$(OPT_VER).o \
		find_class_borders$(OPT_VER).o \
		agf_two_class$(OPT_VER).o \
		agf_util$(OPT_VER).o \
		agf_eval$(OPT_VER).o \
		agf_calibrate$(OPT_VER).o \
		agf_io$(OPT_VER).o \
		agf_fconv$(OPT_VER).o \
		knn_lib$(OPT_VER).o \
		cluster_lib$(OPT_VER).o \
		multi_parse$(OPT_VER).o \
		multiclass_hier$(OPT_VER).o \
		multiclass$(OPT_VER).o \
		multiclass_c$(OPT_VER).o \
		multi_control$(OPT_VER).o \
		multi_control_experimental$(OPT_VER).o \
		multiclass_methods$(OPT_VER).o \
		classifier_obj$(OPT_VER).o \
		binaryclassifier$(OPT_VER).o \
		agf2class$(OPT_VER).o \
		direct_classifier$(OPT_VER).o \
		svmkernel$(OPT_VER).o \
		svm2class$(OPT_VER).o \
		linearclassifier$(OPT_VER).o \
		svm_multi$(OPT_VER).o \
		batch_borders$(OPT_VER).o \
		solve_cond_prob$(OPT_VER).o

INCLUDE_FILES = agf_util.h \
		agf_eval.h \
		agf_calibrate.h \
		agf_metric2.h \
		agf_engine.h \
		knn_lib.h \
		cluster_lib.h \
		agf_io.h \
		agf_fconv.h \
		agf_defs.h \
		agf_two_class.h \
		sample_class_borders.h \
		find_class_borders.h \
		agf_calc_pdf.h \
		adgaf.h \
		classifier_obj.h \
		binaryclassifier.h \
		agf2class.h \
		svmkernel.h \
		svm2class.h \
		linearclassifier.h \
		svm_multi.h \
		direct_classifier.h \
		multi_parse.h \
		multiclass.h \
		multiclass_c.h \
		multiclass_hier.h \
		multi_control.h \
		multi_control_experimental.h \
		multiclass_methods.h \
		agf_classify.h \
		batch_borders.h \
		solve_cond_prob.h

LIBRARY_INCLUDE = agf_lib.h

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

test: $(BIN_PATH)/agf_testsuite$(EXE_EXT)
	agf_testsuite

#old code that still might be useful:
OLD2=$(addsuffix $(OPT_VER)$(EXE_EXT), $(OLD))
OLD3=$(addprefix $(BIN_PATH)/, $(addsuffix $(EXE_EXT), $(OLD)))

old: $(OLD2)

install_old: $(OLD3)

$(EXECUTABLES3) $(OLD3): $(BIN_PATH)/%$(EXE_EXT): %$(OPT_VER)$(EXE_EXT)
	install $< $@

clean:
	rm -f $(OBJECT_FILES) 
	rm -f $(EXECUTABLES2) 
	rm -f $(LIBNAME)
	rm -f $(OLD2)

$(EXECUTABLES2) $(OLD2): %$(OPT_VER)$(EXE_EXT): %.cc $(LIBNAME) $(LIBPETEY) $(LIBRARY_INCLUDE)
	$(CC) $(CFLAGS) -o $@ $< $(LIBNAME) $(LD_OPTIONS)

$(LIBNAME): $(OBJECT_FILES) makefile
	ar -rs $(LIBNAME) $(OBJECT_FILES)

$(OBJECT_FILES): %$(OPT_VER).o: %.cc $(LIBRARY_INCLUDE)
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBRARY_INCLUDE): $(INCLUDE_FILES)
	touch $(LIBRARY_INCLUDE)

svm-binary: svm-binary.cc $(LIBSVM)/svm.o $(LIBSVM)/svm.h $(LIBNAME) $(LIBPETEY) $(LIBRARY_INCLUDE)
	$(CC) $(CFLAGS) -I$(LIBSVM) svm-binary.cc -o svm-binary $(LIBSVM)/svm.o

#$(LIBPETEY):
#	make -C $(DATASET_PATH) 

