SHELL=/bin/bash

TIME=/bin/time

#minimum number of border samples:
MINSAMPLE=5
#maximum number of border samples:
MAXSAMPLE=200
NTEST=20
NTRIAL=5
#fraction of data used for testing:
FRAC=0.4

#number of trials for collecting stats:
NRUN=10

#which dataset to use:
#VER=heart
#VER=shuttle
VER=covtype

#number of classes in each dataset:
NCLS_heart=2
NCLS_shuttle=7
NCLS_covtype=7
NCLS_segment=7
NCLS_sat=6
NCLS_vehicle=4
NCLS_poker=10
NCLS_mnist=10
NCLS_letter=26
NCLS_codrna=2
NCLS_protein=3
NCLS_pendigits=10
NCLS_dna=3
NCLS_usps=10
NCLS_seismic=3
NCLS_ijcnn1=2
NCLS_leu=2
NCLS_mushrooms=2
NCLS_madelon=2
NCLS_phishing=2
NCLS_splice=2
#this one actually sets the number of classes:
NCLS_humidity=8

NCLS=$(NCLS_$(VER))

ifeq ($(NCLS),2)
  AGF_TRAIN=class_borders
  BORDERS_TEST=classify_b
  ACCELERATE_COMMAND=class_borders -Z
  ACCELERATOR2=$(ACCELERATOR).brd
  AGFMODEL2=$(AGFMODEL).brd
else
  AGF_TRAIN=multi_borders
  AGF_PREFIX=print_control -Q 5 $(NCLS) |
  BORDERS_TEST=classify_m
  ACCELERATE_COMMAND=svm_accelerate
  AGFBASE=$(AGFMODEL)
  ACCTRAIN=$(TRAIN)
  ACCELERATOR2=$(ACCELERATOR)
  AGFMODEL2=$(AGFMODEL)
endif

#special pre-processing:
#PRE_mnist=-nS 20

ifdef PRE_$(VER)
  PRE=$(PRE_$(VER))
else
  PRE=-n
endif

#options for svm-train:
SVMOPT_heart=-c 0.5 -g 0.01
SVMOPT_sat=-c 50 -g 0.1
SVMOPT_segment=-c 100 -g 0.1 
SVMOPT_vehicle=-c 100 -g 0.01
SVMOPT_pendigits=-c 50 -g 0.01
SVMOPT_letter=-c 50 -h 0
SVMOPT_seismic=-c 50 -h 0
SVMOPT_humidity=-h 0 -c 50
SVMOPT_mnist=-c 50
SVMOPT_mushrooms=-c 50 -h 0
SVMOPT_codrna=-h 0

SVMOPT=$(SVMOPT_$(VER))

#options for lin-train:
LINOPT_sat=-s 4

LINOPT=$(LINOPT_$(VER))

#number of nearest neighbours for knn:
K_heart=-k 11
K_shuttle=-k 1
K_covtype=-k 11
K_segment=-k 3
K_sat=-k 5
K_vehicle=-k 9
K_pendigits=-k 1
K_usps=-k 1
K_letter=-k 1
K_seismic=-k 31
K_phishing=-k 3
K_splice=-k 1
K_humidity=-k 51
K_mnist=-k 9
K_ijcnn1=-k 11
K_mushrooms=-k 3

K=$(K_$(VER))

#number of border samples (or other options) for "accelerated" SVM model:
NSAMP_heart=-s 100
NSAMP_segment=-s 50
NSAMP_vehicle=-s 50
NSAMP_sat=-s 200
NSAMP_letter=-s 75
#NSAMP_codrna=-s 500 -r-0.1
NSAMP_codrna=-s 500
NSAMP_pendigits=-s 50
NSAMP_dna=-s 500
NSAMP_usps=-s 50
NSAMP_seismic=-s 200
NSAMP_ijcnn1=-s 500
NSAMP_phishing=-s 500
#NSAMP_splice=-s 200 -r-0.99
#NSAMP_splice=-s 1000
NSAMP_splice=-s 500
NSAMP_humidity=-s 200
NSAMP_mnist=-s 500
NSAMP_mushrooms=-s 200

NSAMP=$(NSAMP_$(VER))

#clean targets for collecting stats: 
#(which parts of the makefile do we want to rerun?)
CLEAN_letter=clean_but_borders
CLEAN_ijcnn1=clean_borders
CLEAN_seismic=clean_borders
CLEAN_codrna=clean_borders
CLEAN_humidity=clean_borders
CLEAN_mnist=clean_borders
ifdef CLEAN_$(VER)
  CLEAN_TARGET=$(CLEAN_$(VER))
else
  #since most datasets are already separated into "test" and "training":
  CLEAN_TARGET=clean
endif

#subsampling:
ifdef SUBSAMPLE_SV
  BUILD_TARGET=build_svm
  CLEAN_TARGET=delete_svsub
  STATFILE=$(VER).sub.txt
else
  #if you make stats results are stuck in here:
  STATFILE=$(VER).txt
endif

ifdef SUBSAMPLE
  BUILD_TARGET=build_svm
  CLEAN_TARGET=clean
  STATFILE=$(VER).sub.txt
else
  #if you make stats results are stuck in here:
  STATFILE=$(VER).txt
endif

#for all datasets except shuttle we keep class ratios constant:
ifneq ($(VER),shuttle)
  CONSTCLASS=C
endif

#where to find scripts (for testing acc. vs. border samples):
SCRIPTDIR=../../scripts

#directory containing data:
DATADIR=../data/$(VER)

#working directory for processing data:
WORKDIR=/cygdrive/e/tmpdata/accel_revision_work/$(VER)

#file names:
BASE=$(WORKDIR)/$(VER)
TRAIN=$(BASE).trn
TEST=$(BASE).tst

#append all results in here:
RESULTS=$(WORKDIR)/results.txt

#one-vs-one control file:
MULTICONTROL=$(WORKDIR)/onevsone.mbc

#model files: SVM, AGF, SVM-borders
ifdef SUBSAMPLE_SV
  SVMMODEL=$(BASE).svmmod.1
  SVMMODEL0=$(BASE).svmmod
else
  SVMMODEL=$(BASE).svmmod
  SVMMODEL0=$(SVMMODEL)
endif

LINMODEL=$(BASE).lin
ACCELERATOR=$(BASE).accel

NBORDER_AGF=$(WORKDIR)/nborder_agf.txt
NBORDER_SVM=$(WORKDIR)/nborder_svm.txt
NSUPPORT=$(WORKDIR)/nsupport.txt

#need this to normalize test data same as training data, if applicable:
ifneq ($(strip $(PRE)),)
  NORMFILE=$(TRAIN).std
  NORMCLAUSE=-a $(NORMFILE)
endif

#phony targets:
.PHONY: all
.PHONY: clean
.PHONY: nborders
.PHONY: stats
.PHONY: subsample_sv

#output results:
all: $(RESULTS)
	@echo "KNN:"
	@echo
	@cat $(WORKDIR)/knn_time.txt
	@cls_comp_stats $(TEST).cls $(TEST).knn
	@echo
	@echo "Linear:"
	@echo
	@cat $(WORKDIR)/lin_ttime.txt
	@cat $(WORKDIR)/lin_time.txt
	@echo
	@cls_comp_stats $(TEST).cls $(TEST).lin
	@echo
	@echo "SVM:"
	@echo
	@cat $(WORKDIR)/svm_ttime.txt
	@cat $(WORKDIR)/svm_time.txt
	@echo
	@cls_comp_stats $(TEST).cls $(TEST).svm
	@echo
	@echo "SVM-borders:"
	@echo
	@cat $(WORKDIR)/accel_ttime.txt
	@cat $(WORKDIR)/acc_time.txt
	@echo
	@cls_comp_stats $(TEST).cls $(TEST).acc
	@echo
	@echo "Comparison:"
	@echo
	@cls_comp_stats -q 0 $(TEST).svm.cls $(TEST).acc
	@echo
	@agf_correlate $(TEST).svm.con $(TEST).acc.con
	@echo

#how many lines does results file contain so far?
NLINE!=cat $(STATFILE) | wc -l
stats: ../calc_stats
	#for I in {$(NLINE)..$(NRUN)}; do 
	for ((I=$(NLINE); I<$(NRUN); I++)); do \
		make $(CLEAN_TARGET); \
		make RESULTS=$(STATFILE) $(BUILD_TARGET); \
	done
	./calc_stats < $(STATFILE)

$(RESULTS): $(TEST).knn.cls $(TEST).lin.cls $(TEST).svm.cls $(TEST).acc.cls
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/knn_time.txt | sed '1q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/knn_time.txt | sed '2q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/lin_ttime.txt | sed '1q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/lin_ttime.txt | sed '2q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/lin_time.txt | sed '1q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/lin_time.txt | sed '2q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/svm_ttime.txt | sed '1q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/svm_ttime.txt | sed '2q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/svm_time.txt | sed '1q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/svm_time.txt | sed '2q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/accel_ttime.txt | sed '1q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/accel_ttime.txt | sed '2q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/acc_time.txt | sed '1q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/acc_time.txt | sed '2q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	cls_comp_stats -Hb $(TEST).cls $(TEST).knn | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	cls_comp_stats -Hb $(TEST).cls $(TEST).lin | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	cls_comp_stats -Hb $(TEST).cls $(TEST).svm | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	cls_comp_stats -Hb $(TEST).cls $(TEST).acc | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	cls_comp_stats -Hb $(TEST).svm.cls $(TEST).acc | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	agf_correlate $(TEST).svm.con $(TEST).acc.con | awk '{print $$3}' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep total_sv $(SVMMODEL) | awk '{print $$2}' >> $(RESULTS)


nborders: $(NBORDER_SVM)

clean:
	rm -f $(WORKDIR)/*

#clean only AGF and SVM-borders:
clean_borders:
	rm -f $(AGFMODEL)*
	rm -f $(ACCELERATOR)*

#clean everything but AGF-borders (get stats on timing at least...):
clean_but_borders:
	rm -f $(ACCELERATOR)
	rm -f $(SVMMODEL)
	rm -f $(TEST).knn.*

SVOLD=$(SVMMODEL).old

#new method of sub-sampling:
delete_svsub:
	rm -f $(SVMMODEL)

ifdef SUBSAMPLE_SV
$(SVMMODEL): $(SVMMODEL0)
	subsample_sv -f $(SUBSAMPLE_SV) $(SVMMODEL0) $(SVMMODEL)
endif

#build SVM only (for sub-sampling):
build_svm: $(TEST).svm.cls $(TEST).svm.con
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/svm_ttime.txt | sed '1q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/svm_ttime.txt | sed '2q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/svm_time.txt | sed '1q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/svm_time.txt | sed '2q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	cls_comp_stats -Hb $(TEST).cls $(TEST).svm | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	wc -l $(TRAIN).svm | awk '{print $$1}' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep total_sv $(SVMMODEL) | awk '{print $$2}' >> $(RESULTS)
	@echo "SVM:"
	@echo
	@cat $(WORKDIR)/svm_ttime.txt
	@cat $(WORKDIR)/svm_time.txt
	@echo
	@cls_comp_stats $(TEST).cls $(TEST).svm
	@echo

#testing:
$(TEST).knn.cls $(TEST).knn.con: $(TRAIN).vec $(TRAIN).cls $(TEST).vec
	$(TIME) -o $(WORKDIR)/knn_time.txt knn $K classify $(TRAIN) $(TEST).vec $(TEST).knn

$(TEST).svm.cls $(TEST).svm.con: $(SVMMODEL) $(TEST).svm
	$(TIME) -o $(WORKDIR)/svm_time.txt svm-predict -b 1 $(TEST).svm $(SVMMODEL) $(TEST).svm.txt
	svmout2agf $(TEST).svm.txt $(TEST).svm

$(TEST).lin.cls $(TEST).lin.con: $(LINMODEL) $(TEST).svm
	$(TIME) -o $(WORKDIR)/lin_time.txt lin-predict $(TEST).svm $(LINMODEL) $(TEST).lin.txt
	svmout2agf $(TEST).lin.txt $(TEST).lin

$(TEST).acc.cls $(TEST).acc.con: $(ACCELERATOR2) $(TEST).vec
	$(TIME) -o $(WORKDIR)/acc_time.txt $(BORDERS_TEST) $(ACCELERATOR) $(TEST).vec $(TEST).acc > $(TEST).acc.txt

#find number of border samples SVM-borders:
$(NBORDER_SVM): $(SVMMODEL) $(TRAIN).vec $(TRAIN).cls
	$(SCRIPTDIR)/mb_nb.sh -g -q $(NTEST) -N $(NTRIAL) -M \
		-s $(MINSAMPLE) -S $(MAXSAMPLE) \
		$(SVMMODEL) $(TRAIN) $(NBORDER_SVM)

#find number of support vectors:
$(NSUPPORT): $(TRAIN).svm
	$(SCRIPTDIR)/svm_nsv.sh -G -q $(NTEST) $(SVMOPT) $(TRAIN) $(NSUPPORT)

#svm-borders "accelerated" svm model
$(ACCELERATOR2): $(SVMMODEL0) $(TRAIN).vec $(TRAIN).cls
	$(TIME) -o $(WORKDIR)/accel_ttime.txt $(ACCELERATE_COMMAND) $(NSAMP) $(SVMMODEL0) $(TRAIN) $(ACCELERATOR)

#svm model:
$(SVMMODEL0): $(TRAIN).svm
	$(TIME) -o $(WORKDIR)/svm_ttime.txt svm-train -b 1 $(SVMOPT) $(TRAIN).svm $(SVMMODEL0)

#agf borders model:
$(LINMODEL): $(TRAIN).svm
	$(TIME) -o $(WORKDIR)/lin_ttime.txt lin-train $(LINOPT) $(TRAIN).svm $(LINMODEL)

#only used for one script:
$(MULTICONTROL):
	print_control -Q 5 $(NCLS) > $(MULTICONTROL)

#make the working directory:
$(WORKDIR)/$(VER).timestamp:
	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	date > $(WORKDIR)/$(VER).timestamp

#convert to LIBSVM format:
$(TRAIN).svm: $(TRAIN).cls $(TRAIN).vec
	agf2ascii -M $(TRAIN) $(TRAIN).svm

$(TEST).svm: $(TEST).cls $(TEST).vec
	agf2ascii -M $(TEST) $(TEST).svm

$(TEST).vec $(TEST).cls: ../work/$(VER)/$(VER).tst.vec $(WORKDIR)/$(VER).timestamp
	make -C .. VER=$(VER) work/$(VER)/$(VER).tst.vec
	cp ../work/$(VER)/$(VER).tst.vec $(TEST).vec
	cp ../work/$(VER)/$(VER).tst.cls $(TEST).cls

$(TRAIN).vec $(TRAIN).cls: ../work/$(VER)/$(VER).trn.vec $(WORKDIR)/$(VER).timestamp
	make -C .. VER=$(VER) work/$(VER)/$(VER).trn.vec
	cp ../work/$(VER)/$(VER).trn.vec $(TRAIN).vec
	cp ../work/$(VER)/$(VER).trn.cls $(TRAIN).cls

