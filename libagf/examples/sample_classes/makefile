#---------------------------------------------------------------------#
#            The user should modify the following macros              #
#            to suit his or her system:                               #
#---------------------------------------------------------------------#

BASE_PATH = /home/lenovo

INCLUDE_PATH = $(BASE_PATH)/include
LIB_PATH =$(BASE_PATH)/lib
BIN_PATH =$(BASE_PATH)/bin

#GSL_INCLUDE=$(INCLUDE_PATH)
GSL_INCLUDE=/usr/local/include
#GSL_LIB=$(LIB_PATH)
GSL_LIB=/usr/local/lib

CC = g++
OPT_VER = -g
CFLAGS = $(OPT_VER) -I$(GSL_INCLUDE) -I$(INCLUDE_PATH)

LD_FLAGS = -L$(GSL_LIB) -L$(LIB_PATH) -lagf$(OPT_VER) -lpetey$(OPT_VER) -lgsl -lgslcblas

#---------------------------------------------------------------------#
#                      End modifiable macros                          #
#                                                                     #
#---------------------------------------------------------------------#

#VPATH=$(LIB_PATH):$(INCLUDE_PATH):..

LIBNAME = libagf$(OPT_VER).a

SAMPLE_CLASS = sample_class classify_sc pdf_sc1 pdf_sc2 sc_borders class2_errvh

OBJECT_FILES = sample_class1_obj.o sample_class2_obj.o

all: $(SAMPLE_CLASS)

install:
	install $(FILE_UT) $(BIN_PATH)
	install $(SAMPLE_CLASS) $(BIN_PATH)

#think I want to move these to a separate directory:
file_conversion: $(FILE_UT)

clean_all: clean comp_clean pdf_clean clean_test

$(SAMPLE_CLASS): % : %.cc $(OBJECT_FILES)
	$(CC) $(CFLAGS) -o $@ $< $(OBJECT_FILES) -lagf$(OPT_VER) $(LD_FLAGS)

$(OBJECT_FILES): %.o : %.cc %.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o
	rm -f $(SAMPLE_CLASS)
	rm -f $(FILE_UT)

########  This part of the makefile is dedicated to testing ########

test: sctest.cls sctest.aa.cls sctest.aa.con sctest.a.cls \
		sctest.a.con sctest.b.cls sctest.b.con
	@echo
	@echo "Truth versus analytic"
	@echo
	cls_comp_stats sctest.cls sctest.aa
	@echo
	@echo "Truth versus AGF"
	@echo
	cls_comp_stats sctest.cls sctest.a
	@echo
	@echo "Truth versus AGF borders"
	@echo
	cls_comp_stats sctest.cls sctest.b
	@echo
	@echo "Analytic versus AGF"
	@echo
	cls_comp_stats sctest.aa.cls sctest.a
	@echo
	@echo "AGF versus AGF borders"
	@echo
	cls_comp_stats sctest.a.cls sctest.b

clean_test:
	rm -f sctest.vec sctest.cls sctest.std
	rm -f sctest.*.con sctest.*.cls sctrain.aa.txt
	rm -f sctrain.vec sctrain.cls sctrain.brd sctrain.bgd sctrain.std
	rm -f sctrain.txt
	rm -f sctest*.txt
	rm -f time_*.txt
	rm -f sctest.*.r

sctest.vec sctest.cls sctest.txt: sample_class
	./sample_class -R 1000 2000 sctest > sctest.txt

sctest.aa.cls sctest.aa.con sctest.aa.txt: sctest.vec classify_sc
	time -o time_analytic.txt \
		./classify_sc sctest.vec sctest.aa 
	C2R sctest.aa sctest.aa.r

sctest.a.cls sctest.a.con time_classify_a.txt: sctest.vec sctrain.vec sctrain.cls
	time -o time_classify_a.txt \
		agf -k 1000 -W 100 classify sctrain sctest.vec sctest.a 
	C2R sctest.a sctest.a.r

sctest.b.cls sctest.b.con time_classify_b.txt: sctest.vec sctrain.brd sctrain.bgd #C2R
	time -o time_classify_b.txt \
		classify_b sctrain sctest.vec sctest.b
	C2R sctest.b sctest.b.r

sctrain.brd sctrain.bgd time_train_b.txt: sctrain.vec sctrain.cls
	time -o time_train_b.txt \
		class_borders -r -0.1 -k 1000 -W 100 -s 250 sctrain sctrain

sctrain.vec sctrain.cls sctrain.txt: sample_class
	./sample_class 5000 10000 sctrain > sctrain.txt 

####  This part of the makefile is dedicated to a comparison with Kohonen's methods  #####

compare:	sctest.cls sctest.aa.con \
		time_knn.txt \
		sctest.knn.cls sctest.knn.con \
		time_classify_a.txt time_classify_b.txt time_train_b.txt \
		sctest.a.cls sctest.a.con sctest.b.cls sctest.b.con \
		time_lvq_train1.txt time_lvq_train2.txt time_lvq.txt \
		sctest.lvq.cls \
		time_svm_train.txt time_svm.txt \
		sctest.svm.cls sctest.svm.con

	@echo KNN performance
	@echo __________________________________
	@echo time:
	@more time_knn.txt
	@echo -e "\nclassification accuracy:"
	@cls_comp_stats sctest.cls sctest.knn
	@echo -e "\naccuracy of conditional probabilies:"
	@agf_correlate sctest.aa.r sctest.knn.r
	@echo

	@echo AGF performance
	@echo __________________________________
	@echo time:
	@more time_classify_a.txt
	@echo -e "\nclassification accuracy:"
	@cls_comp_stats sctest.cls sctest.a
	@echo -e "\naccuracy of conditional probabilies:"
	@agf_correlate sctest.aa.r sctest.a.r
	@echo

	@echo AGF borders performance
	@echo __________________________________
	@echo training time:
	@more time_train_b.txt
	@echo -e "\nclassification time"
	@more time_classify_b.txt
	@echo -e "\nclassification accuracy:"
	@cls_comp_stats sctest.cls sctest.b
	@echo -e "\naccuracy of conditional probabilies:"
	@agf_correlate sctest.aa.r sctest.b.r
	@echo

	@echo LVQ performance
	@echo __________________________________
	@echo training time:
	@more time_lvq_train1.txt
	@more time_lvq_train2.txt
	@echo -e "\nclassification time:"
	@more time_lvq.txt
	@echo -e "\nclassification accuracy:"
	@cls_comp_stats sctest.cls sctest.lvq
	@echo

	@echo SVM performance
	@echo __________________________________
	@echo training time:
	@more time_svm_train.txt
	@echo -e "\nclassification time"
	@more time_svm.txt
	@echo -e "\nclassification accuracy:"
	@cls_comp_stats sctest.cls sctest.svm
	@echo -e "\naccuracy of conditional probabilies:"
	@agf_correlate sctest.aa.r sctest.svm.r
	@echo

comp_clean:
	rm -f sctest.knn.*
	rm -f sctrain.lvq.cod
	rm -f sctest.lvq.*
	#rm -f time_*.txt
	rm -f sctest.svm*
	rm -f sctrain.svm*

# K-nearest neighbours method:
sctest.knn.cls sctest.knn.con time_knn.txt: sctrain.vec sctrain.cls sctest.vec #C2R
	time -o time_knn.txt knn -k 21 classify sctrain sctest.vec sctest.knn
	C2R sctest.knn sctest.knn.r

# LVQ training:
sctrain.lvq.cod time_lvq_train1.txt time_lvq_train2.txt: sctrain.txt
	time -o time_lvq_train1.txt \
		propinit -din sctrain.txt -cout sctrain.lvq.cod -noc 1000 -knn 5
	time -o time_lvq_train2.txt \
		olvq1 -din sctrain.txt -cin sctrain.lvq.cod -cout sctrain.lvq.cod \
			-rlen 75000 -alpha 0.1

# LVQ classification:
sctest.lvq.cls time_lvq.txt: sctrain.lvq.cod sctest.txt #lvq2bin
	time -o time_lvq.txt \
		classify -din sctest.txt -cin sctrain.lvq.cod -dout sctest.lvq.txt 
	lvq2agf sctest.lvq.txt sctest.lvq

# SOM training:
sctest.som.cod time_som_train.txt: sctrain.txt sc_som_train.bash
	time -o time_som_train.txt sc_som_train.bash

# SOM classification:

# SVM training:
sctrain.svm.model time_svm_train.txt: sctrain.vec sctrain.cls #bincls2svm
	agf2ascii -M sctrain > sctrain.svm
	time -o time_svm_train.txt svm-train -t 2 -g 0.5 -s 0 -b 1 -c 100 sctrain.svm

#SVM classification:
sctest.svm.cls sctest.svm.con time_svm.txt: sctrain.svm.model sctest.vec \
		sctest.cls #bincls2svm svmout2agf C2R
	agf2ascii -M sctest > sctest.svm
	time -o time_svm.txt svm-predict -b 1 sctest.svm sctrain.svm.model \
		sctest.svmout
	svmout2agf sctest.svmout sctest.svm
	C2R sctest.svm sctest.svm.r


########  This part of the makefile is dedicated to testing pdf calculation methods ########

test_pdf:  cls1_pdf.dat cls2_pdf.dat cls1_knn_pdf.dat cls1_agf_pdf.dat \
		cls2_knn_pdf.dat cls2_agf_pdf.dat
	@echo Class 1:
	@echo
	@echo Nearest-neighbours accuracy:
	@agf_correlate cls1_pdf.dat cls1_knn_pdf.dat
	@echo 
	@echo AGF accuracy:
	@agf_correlate cls1_pdf.dat cls1_agf_pdf.dat
	@echo 
	@echo Class 2:
	@echo
	@echo Nearest-neighbours accuracy:
	@agf_correlate cls2_pdf.dat cls2_knn_pdf.dat
	@echo 
	@echo AGF accuracy:
	@agf_correlate cls2_pdf.dat cls2_agf_pdf.dat

pdf_clean:
	rm -f cls1_*.dat
	rm -f cls2_*.dat
	rm -f cls?_train.cls
	rm -f cls?_train.vec
	rm -f cls?_train.txt
	rm -f cls1_*pdf.txt
	rm -f cls2_*pdf.txt

cls1_pdf.dat: sctest.vec pdf_sc1
	./pdf_sc1 sctest.vec cls1_pdf.dat > cls1_pdf.txt

cls2_pdf.dat: sctest.vec pdf_sc2
	./pdf_sc2 sctest.vec cls2_pdf.dat > cls2_pdf.txt

cls1_knn_pdf.dat: sctest.vec cls1_train.vec
	time knn -k 50 pdf cls1_train sctest.vec cls1_knn_pdf > cls1_knn_pdf.txt

cls1_agf_pdf.dat: sctest.vec cls1_train.vec
	time agf -V 1 -k 1000 -W 100 pdf cls1_train sctest.vec cls1_agf_pdf > cls1_agf_pdf.txt

cls2_knn_pdf.dat: sctest.vec cls2_train.vec
	time knn -k 100 pdf cls2_train sctest.vec cls2_knn_pdf > cls2_knn_pdf.txt

cls2_agf_pdf.dat: sctest.vec cls2_train.vec
	time agf -k 1000 -W 100 pdf cls2_train sctest.vec cls2_agf_pdf > cls2_adg_pdf.txt

cls1_train.vec: sample_class
	./sample_class 5000 0 cls1_train > cls1_train.txt

cls2_train.vec: sample_class
	./sample_class 0 10000 cls2_train > cls1_train.txt

