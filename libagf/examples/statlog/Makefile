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
  AGF_PREFIX=print_control -Q 6 $(NCLS) |
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

#options for AGF-borders:
AGFOPT_heart=-W 10
AGFOPT_shuttle=-W 1.5 -k 200 -s 250
AGFOPT_covtype=
AGFOPT_segment=-W 3 -s 100 -k 150
AGFOPT_sat=-W 5 -s 200 -k 100
AGFOPT_poker=-k 200 -W 2
AGFOPT_mnist=-k 400 -W 25
AGFOPT_vehicle=-W 5 -s 100 -k 100
AGFOPT_letter=-W 3 -k 200 -s 1000
#AGFOPT_codrna=-k 200 -W 10 -s 500 -r-0.28
AGFOPT_codrna=-k 200 -W 10 -s 500
AGFOPT_protein=-k 200 -W 10
AGFOPT_pendigits=-k 200 -W 2 -s 200
AGFOPT_dna=-k 300 -s 500 -W 20
AGFOPT_usps=-k 100 -W 1 -s 250
AGFOPT_seismic=-W 40 -k 500 -s 200
AGFOPT_ijcnn1=-W 10 -k 400 -s 500
AGFOPT_madelon=-W 10 -k 200
AGFOPT_phishing=-W 7 -k 250 -s 500
#AGFOPT_splice=-W 15 -s 200 -r-0.23
AGFOPT_splice=-W 15 -s 500
#AGFOPT_splice=-W 15 -s 1000
AGFOPT_humidity=-k 400 -W 10 -s 200
AGFOPT_mnist=-k 200 -W 1 -s 500
AGFOPT_mushrooms=-k 300 -W 3 -s 500

AGFOPT=$(AGFOPT_$(VER))

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
DATADIR=data/$(VER)

#working directory for processing data:
WORKDIR=work/$(VER)

#file names:
BASE=$(WORKDIR)/$(VER)
TRAIN=$(BASE).trn
TEST=$(BASE).tst

#append all results in here:
RESULTS=$(WORKDIR)/results.txt

#one-vs-one control file:
MULTICONTROL=$(WORKDIR)/onevsone.mbc

#model files: SVM, AGF, SVM-borders
SVMMODEL=$(BASE).svmmod
AGFMODEL=$(BASE).agf
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

#output results:
all: $(RESULTS)
	@echo "KNN:"
	@echo
	@cat $(WORKDIR)/knn_time.txt
	@cls_comp_stats $(TEST).cls $(TEST).knn
	@echo
	@echo "AGF-borders:"
	@echo
	@cat $(WORKDIR)/agf_ttime.txt
	@cat $(WORKDIR)/agf_time.txt
	@echo
	@cls_comp_stats $(TEST).cls $(TEST).agf
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
stats: calc_stats
	#for I in {$(NLINE)..$(NRUN)}; do 
	for ((I=$(NLINE); I<$(NRUN); I++)); do \
		make $(CLEAN_TARGET); \
		make RESULTS=$(STATFILE) $(BUILD_TARGET); \
	done
	./calc_stats < $(STATFILE)

$(RESULTS): $(TEST).knn.cls $(TEST).agf.cls $(TEST).svm.cls $(TEST).acc.cls
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/knn_time.txt | sed '1q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/knn_time.txt | sed '2q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/agf_ttime.txt | sed '1q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/agf_ttime.txt | sed '2q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/agf_time.txt | sed '1q;d' | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	grep -o "[0-9]*\.[0-9]*" $(WORKDIR)/agf_time.txt | sed '2q;d' | tr -d '\n' >> $(RESULTS)
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
	cls_comp_stats -Hb $(TEST).cls $(TEST).agf | tr -d '\n' >> $(RESULTS)
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

$(VEHICLE_TEST):
	make VER=vehicle VEHICLE_TEST=$(VEHICLE_TEST)
	mv work/vehicle/vehicle.tst.knn.cls work/vehicle/vehicle.$(VEHICLE_TEST).tst.knn.cls
	mv work/vehicle/vehicle.tst.knn.con work/vehicle/vehicle.$(VEHICLE_TEST).tst.knn.con
	mv work/vehicle/knn_time.txt work/vehicle/knn_time.$(VEHICLE_TEST).txt
	mv work/vehicle/vehicle.tst.agf.cls work/vehicle/vehicle.$(VEHICLE_TEST).tst.agf.cls
	mv work/vehicle/vehicle.tst.agf.con work/vehicle/vehicle.$(VEHICLE_TEST).tst.agf.con
	mv work/vehicle/agf_ttime.txt work/vehicle/agf_ttime.$(VEHICLE_TEST).txt
	mv work/vehicle/agf_time.txt work/vehicle/agf_time.$(VEHICLE_TEST).txt
	mv work/vehicle/vehicle.tst.svm.cls work/vehicle/vehicle.$(VEHICLE_TEST).tst.svm.cls
	mv work/vehicle/vehicle.tst.svm.con work/vehicle/vehicle.$(VEHICLE_TEST).tst.svm.con
	mv work/vehicle/svm_ttime.txt work/vehicle/svm_ttime.$(VEHICLE_TEST).txt
	mv work/vehicle/svm_time.txt work/vehicle/svm_time.$(VEHICLE_TEST).txt
	mv work/vehicle/vehicle.tst.acc.cls work/vehicle/vehicle.$(VEHICLE_TEST).tst.acc.cls
	mv work/vehicle/vehicle.tst.acc.con work/vehicle/vehicle.$(VEHICLE_TEST).tst.acc.con
	mv work/vehicle/accel_ttime.txt work/vehicle/accel_ttime.$(VEHICLE_TEST).txt
	mv work/vehicle/acc_time.txt work/vehicle/acc_time.$(VEHICLE_TEST).txt
	rm $(TEST).cls $(TEST).vec
	rm $(TRAIN).cls $(TRAIN).vec

#n-fold cross validation for vehicle dataset:
#(at some point we hope to make this more elegant and to work with all the other
#datasets that don't have separate training sets)
vehicle:
	make xaa VEHICLE_TEST=xaa VER=vehicle
	make xab VEHICLE_TEST=xab VER=vehicle
	make xac VEHICLE_TEST=xac VER=vehicle
	make xad VEHICLE_TEST=xad VER=vehicle
	make xae VEHICLE_TEST=xae VER=vehicle
	make xaf VEHICLE_TEST=xaf VER=vehicle
	make xag VEHICLE_TEST=xag VER=vehicle
	make xah VEHICLE_TEST=xah VER=vehicle
	make xai VEHICLE_TEST=xai VER=vehicle
	cat work/vehicle/vehicle.xa?.tst.knn.cls > work/vehicle/vehicle.all.tst.knn.cls
	cat work/vehicle/vehicle.xa?.tst.knn.con > work/vehicle/vehicle.all.tst.knn.con
	cat work/vehicle/vehicle.xa?.tst.agf.cls > work/vehicle/vehicle.all.tst.agf.cls
	cat work/vehicle/vehicle.xa?.tst.agf.con > work/vehicle/vehicle.all.tst.agf.con
	cat work/vehicle/vehicle.xa?.tst.svm.cls > work/vehicle/vehicle.all.tst.svm.cls
	cat work/vehicle/vehicle.xa?.tst.svm.con > work/vehicle/vehicle.all.tst.svm.con
	cat work/vehicle/vehicle.xa?.tst.acc.cls > work/vehicle/vehicle.all.tst.acc.cls
	cat work/vehicle/vehicle.xa?.tst.acc.con > work/vehicle/vehicle.all.tst.acc.con
	cat $(addsuffix $(VEHICLE_SUFFIX),$(addprefix data/vehicle/, $(VEHICLE_FILES))) | sed 's/opel/0/;s/saab/1/;s/bus/2/;s/van/3/' | lvq2agf -H work/vehicle/vehicle.all.tst

$(TEST).$(VEHICLE_TEST).%.cls: $(TEST).%.cls
	cp $< $@

$(TEST).$(VEHICLE_TEST).%.vec: $(TEST).%.vec
	cp $< $@

#testing:
$(TEST).knn.cls $(TEST).knn.con: $(TRAIN).vec $(TRAIN).cls $(TEST).vec
	$(TIME) -o $(WORKDIR)/knn_time.txt knn $K classify $(TRAIN) $(TEST).vec $(TEST).knn

$(TEST).svm.cls $(TEST).svm.con: $(SVMMODEL) $(TEST).svm
	$(TIME) -o $(WORKDIR)/svm_time.txt svm-predict -b 1 $(TEST).svm $(SVMMODEL) $(TEST).svm.txt
	svmout2agf $(TEST).svm.txt $(TEST).svm

$(TEST).agf.cls $(TEST).agf.con: $(AGFMODEL2) $(TEST).vec
	$(TIME) -o $(WORKDIR)/agf_time.txt $(BORDERS_TEST) $(AGFMODEL) $(TEST).vec $(TEST).agf > $(TEST).agf.txt

$(TEST).acc.cls $(TEST).acc.con: $(ACCELERATOR2) $(TEST).vec
	$(TIME) -o $(WORKDIR)/acc_time.txt $(BORDERS_TEST) $(ACCELERATOR) $(TEST).vec $(TEST).acc > $(TEST).acc.txt

#find number of border samples SVM-borders:
$(NBORDER_SVM): $(SVMMODEL) $(TRAIN).vec $(TRAIN).cls
	$(SCRIPTDIR)/mb_nb.sh -g -q $(NTEST) -N $(NTRIAL) -M \
		-s $(MINSAMPLE) -S $(MAXSAMPLE) \
		$(SVMMODEL) $(TRAIN) $(NBORDER_SVM)

#find number of border samples AGF-borders:
$(NBORDER_AGF): $(MULTICONTROL) $(TRAIN).vec $(TRAIN).cls
	$(SCRIPTDIR)/mb_nb.sh -g -q $(NTEST) -N $(NTRIAL) \
		-s $(MINSAMPLE) -S $(MAXSAMPLE) $(AGFOPT) \
		$(MULTICONTROL) $(TRAIN) $(NBORDER_AGF)

#find number of support vectors:
$(NSUPPORT): $(TRAIN).svm
	$(SCRIPTDIR)/svm_nsv.sh -G -q $(NTEST) $(SVMOPT) $(TRAIN) $(NSUPPORT)

#svm-borders "accelerated" svm model
$(ACCELERATOR2): $(SVMMODEL) $(TRAIN).vec $(TRAIN).cls
	$(TIME) -o $(WORKDIR)/accel_ttime.txt $(ACCELERATE_COMMAND) $(NSAMP) $(SVMMODEL) $(TRAIN) $(ACCELERATOR)

#svm model:
$(SVMMODEL): $(TRAIN).svm
	$(TIME) -o $(WORKDIR)/svm_ttime.txt svm-train -b 1 $(SVMOPT) $(TRAIN).svm $(SVMMODEL)

#agf borders model:
$(AGFMODEL2): $(TRAIN).svm
	$(AGF_PREFIX) $(TIME) -o $(WORKDIR)/agf_ttime.txt $(AGF_TRAIN) $(AGFOPT) $(TRAIN) $(AGFBASE) $(AGFMODEL)

#only used for one script:
$(MULTICONTROL):
	print_control -Q 6 $(NCLS) > $(MULTICONTROL)

#convert to LIBSVM format:
$(TRAIN).svm: $(TRAIN).cls $(TRAIN).vec
	agf2ascii -M $(TRAIN) $(TRAIN).svm

$(TEST).svm: $(TEST).cls $(TEST).vec
	agf2ascii -M $(TEST) $(TEST).svm

#file conversion:
#steps:
#1. convert to libAGF format and copy to working directory
#2. split into test and training, if applicable
#3. pre-processing: 
#     - normalization (optional)
#     - SVD (optional)
#     - make class labels go from 0 to Nc-1
#list of datasets from LIBSVM website with test and training:

#list of datasets from LIBSVM website with test and training:
SVMTTDATA=madelon \
	ijcnn1 \
	seismic \
	usps \
	dna \
	pendigits \
	codrna \
	poker 
	#splice

#list of datasets that need to be separated into test and training:
ALLINONE=	heart \
		segment \
		mushrooms \
		phishing \
		covtype \
		letter \
		humidity

#make the working directory:
$(WORKDIR)/$(VER).timestamp:
	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	date > $(WORKDIR)/$(VER).timestamp

#UCI datasets in LVQ format (sans header) with only training data:
$(WORKDIR)/heart.cls: $(WORKDIR)/%.cls: $(DATADIR)/%.dat.txt $(WORKDIR)/$(VER).timestamp
	lvq2agf -H $< $(WORKDIR)/$(VER)

$(WORKDIR)/segment.cls: $(DATADIR)/segment.dat.txt $(WORKDIR)/$(VER).timestamp
	lvq2agf -H $< $(WORKDIR)/$(VER)
	#3rd feature must be removed since it only takes on one value:
	agf_preprocess -F $(WORKDIR)/$(VER) $(WORKDIR)/segment 0 1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18

$(WORKDIR)/letter.cls: $(DATADIR)/letter-recognition.data.txt $(WORKDIR)/$(VER).timestamp
	sed 's/A/0/;s/B/1/;s/C/2/;s/D/3/;s/E/4/;s/F/5/;s/G/6/;s/H/7/;s/I/8/;s/J/9/;s/K/10/;s/L/11/;s/M/12/;s/N/13/;s/O/14/;s/P/15/;s/Q/16/;s/Q/16/;s/R/17/;s/S/18/;s/T/19/;s/U/20/;s/V/21/;s/W/22/;s/X/23/;s/Y/24/;s/Z/25/;' < $< | awk -F ',' '{print $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9, $$10, $$11, $$12, $$13, $$14, $$15, $$16, $$1}' > $(WORKDIR)/letter.txt
	lvq2agf -H $(WORKDIR)/letter.txt $(WORKDIR)/letter

#UCI datasets in LVQ format (sans header) with both training and test data:
$(WORKDIR)/sat.trn.raw.cls $(WORKDIR)/shuttle.trn.raw.cls: $(WORKDIR)/%.trn.raw.cls: $(DATADIR)/%.trn.txt $(WORKDIR)/$(VER).timestamp
	lvq2agf -H $(DATADIR)/$(VER).trn.txt $(TRAIN).raw

$(WORKDIR)/sat.tst.raw.cls $(WORKDIR)/shuttle.tst.raw.cls: $(WORKDIR)/%.tst.raw.cls: $(DATADIR)/%.tst.txt $(WORKDIR)/$(VER).timestamp
	lvq2agf -H $(DATADIR)/$(VER).tst.txt $(TEST).raw

#LIBSVM datasets with only training data:
$(WORKDIR)/mushrooms.cls $(WORKDIR)/phishing.cls $(WORKDIR)/covtype.cls: $(WORKDIR)/%.cls: $(DATADIR)/% $(WORKDIR)/$(VER).timestamp
	svm2agf -E 0 $(DATADIR)/$(VER) $(WORKDIR)/$(VER)

#LIBSVM datasets with both test and training data:
$(addsuffix .trn.raw.cls, $(addprefix $(WORKDIR)/, $(SVMTTDATA))): $(WORKDIR)/%.trn.raw.cls: $(DATADIR)/% $(WORKDIR)/$(VER).timestamp
	svm2agf -E 0 $(DATADIR)/$(VER) $(TRAIN).raw

$(addsuffix .tst.raw.cls, $(addprefix $(WORKDIR)/, $(SVMTTDATA))): $(WORKDIR)/%.tst.raw.cls: $(DATADIR)/%.t $(WORKDIR)/$(VER).timestamp
	svm2agf -E 0 $(DATADIR)/$(VER).t $(TEST).raw

#rearrange splice features:
SED_PROGRAM='s/1\./3\./;s/2\./2\./;s/3\./0\./;s/4\./1\./;'
$(WORKDIR)/splice.trn.raw.cls: $(DATADIR)/splice $(WORKDIR)/$(VER).timestamp
	sed $(SED_PROGRAM) < $< > $(WORKDIR)/splice.txt
	svm2agf -E 0 $(WORKDIR)/splice.txt $(WORKDIR)/splice.trn.raw

$(WORKDIR)/splice.tst.raw.cls: $(DATADIR)/splice.t $(WORKDIR)/$(VER).timestamp
	sed $(SED_PROGRAM) < $< > $(WORKDIR)/splice.t.txt
	svm2agf -E 0 $(WORKDIR)/splice.t.txt $(WORKDIR)/splice.tst.raw

$(WORKDIR)/mnist.trn.raw.cls: $(DATADIR)/mnist $(WORKDIR)/$(VER).timestamp
	svm2agf -E 0 $(DATADIR)/$(VER) $(TRAIN).raw
	agf_preprocess -F $(TRAIN).raw $(TRAIN).raw $(MNIST_INDICES)

$(WORKDIR)/mnist.tst.raw.cls: $(DATADIR)/mnist.t $(WORKDIR)/$(VER).timestamp
	svm2agf -E 0 $(DATADIR)/$(VER).t $(TEST).raw
	#because test data has different dimensions from the training
	#data, feature selection has to be done separately:
	agf_preprocess -F $(TEST).raw $(TEST).raw $(MNIST_INDICES)

#my own dataset:
$(WORKDIR)/humidity.cls: ../humidity_data/train_SQ36s.dat ../humidity_data/est_bt_cld_land_z0.00.vec
	cp ../humidity_data/est_bt_cld_land_z0.00.vec $(WORKDIR)/humidity.vec
	float_to_class -q $(NCLS_humidity) -g ../humidity_data/train_SQ36s.dat $(WORKDIR)/humidity.cls 0.00007 0.001

#convert raw data into normalized training and test sets:

#datasets where all the data is lumped together and needs to be separated:
$(addsuffix .tst.raw.cls, $(addprefix $(WORKDIR)/, $(ALLINONE))): $(WORKDIR)/%.tst.raw.cls: $(WORKDIR)/%.trn.raw.cls

#separate test from training:
$(addsuffix .trn.raw.cls, $(addprefix $(WORKDIR)/, $(ALLINONE))): $(WORKDIR)/%.trn.raw.cls: $(WORKDIR)/%.cls
#$(TRAIN).raw.cls $(TEST).raw.cls: $(WORKDIR)/$(VER).cls
	agf_preprocess $(PRE) -zf $(FRAC) $(WORKDIR)/$(VER) $(TRAIN).raw $(TEST).raw

#process training data:
$(TRAIN).vec $(TRAIN).cls $(NORMFILE): $(TRAIN).raw.cls
	agf_preprocess $(PRE) -U $(TRAIN).raw $(TRAIN)
ifdef SUBSAMPLE
	subsample_special -z$(CONSTCLASS)f $(SUBSAMPLE) $(TRAIN) $(WORKDIR)/dum $(TRAIN)
endif

#process test data:
$(TEST).vec $(TEST).cls: $(TEST).raw.cls $(NORMFILE)
	agf_preprocess $(NORMCLAUSE) -U $(TEST).raw $(TEST)

#we like the mnist dataset because it takes 3 days to train in LIBSVM:
MNIST_INDICES = 33 \
34 \
35 \
36 \
37 \
38 \
39 \
40 \
41 \
42 \
43 \
44 \
45 \
46 \
47 \
48 \
61 \
62 \
63 \
64 \
65 \
66 \
67 \
68 \
69 \
70 \
71 \
72 \
73 \
74 \
75 \
76 \
77 \
78 \
79 \
80 \
88 \
89 \
90 \
91 \
92 \
93 \
94 \
95 \
96 \
97 \
98 \
99 \
100 \
101 \
102 \
103 \
104 \
105 \
106 \
107 \
108 \
109 \
110 \
115 \
116 \
117 \
118 \
119 \
120 \
121 \
122 \
123 \
124 \
125 \
126 \
127 \
128 \
129 \
130 \
131 \
132 \
133 \
134 \
135 \
136 \
137 \
138 \
142 \
143 \
144 \
145 \
146 \
147 \
148 \
149 \
150 \
151 \
152 \
153 \
154 \
155 \
156 \
157 \
158 \
159 \
160 \
161 \
162 \
163 \
164 \
165 \
166 \
170 \
171 \
172 \
173 \
174 \
175 \
176 \
177 \
178 \
179 \
180 \
181 \
182 \
183 \
184 \
185 \
186 \
187 \
188 \
189 \
190 \
191 \
192 \
193 \
194 \
195 \
197 \
198 \
199 \
200 \
201 \
202 \
203 \
204 \
205 \
206 \
207 \
208 \
209 \
210 \
211 \
212 \
213 \
214 \
215 \
216 \
217 \
218 \
219 \
220 \
221 \
222 \
223 \
225 \
226 \
227 \
228 \
229 \
230 \
231 \
232 \
233 \
234 \
235 \
236 \
237 \
238 \
239 \
240 \
241 \
242 \
243 \
244 \
245 \
246 \
247 \
248 \
249 \
250 \
251 \
252 \
253 \
254 \
255 \
256 \
257 \
258 \
259 \
260 \
261 \
262 \
263 \
264 \
265 \
266 \
267 \
268 \
269 \
270 \
271 \
272 \
273 \
274 \
275 \
276 \
277 \
278 \
279 \
280 \
281 \
282 \
283 \
284 \
285 \
286 \
287 \
288 \
289 \
290 \
291 \
292 \
293 \
294 \
295 \
296 \
297 \
298 \
299 \
300 \
301 \
302 \
303 \
304 \
305 \
306 \
307 \
308 \
309 \
310 \
311 \
312 \
313 \
314 \
315 \
316 \
317 \
318 \
319 \
320 \
321 \
322 \
323 \
324 \
325 \
326 \
327 \
328 \
329 \
330 \
331 \
332 \
333 \
334 \
335 \
337 \
338 \
339 \
340 \
341 \
342 \
343 \
344 \
345 \
346 \
347 \
348 \
349 \
350 \
351 \
352 \
353 \
354 \
355 \
356 \
357 \
358 \
359 \
360 \
361 \
362 \
365 \
366 \
367 \
368 \
369 \
370 \
371 \
372 \
373 \
374 \
375 \
376 \
377 \
378 \
379 \
380 \
381 \
382 \
383 \
384 \
385 \
386 \
387 \
388 \
389 \
390 \
394 \
395 \
396 \
397 \
398 \
399 \
400 \
401 \
402 \
403 \
404 \
405 \
406 \
407 \
408 \
409 \
410 \
411 \
412 \
413 \
414 \
415 \
416 \
417 \
418 \
422 \
423 \
424 \
425 \
426 \
427 \
428 \
429 \
430 \
431 \
432 \
433 \
434 \
435 \
436 \
437 \
438 \
439 \
440 \
441 \
442 \
443 \
444 \
445 \
446 \
447 \
449 \
450 \
451 \
452 \
453 \
454 \
455 \
456 \
457 \
458 \
459 \
460 \
461 \
462 \
463 \
464 \
465 \
466 \
467 \
468 \
469 \
470 \
471 \
472 \
473 \
474 \
475 \
478 \
479 \
480 \
481 \
482 \
483 \
484 \
485 \
486 \
487 \
488 \
489 \
490 \
491 \
492 \
493 \
494 \
495 \
496 \
497 \
498 \
499 \
500 \
501 \
502 \
503 \
505 \
506 \
507 \
508 \
509 \
510 \
511 \
512 \
513 \
514 \
515 \
516 \
517 \
518 \
519 \
520 \
521 \
522 \
523 \
524 \
525 \
526 \
527 \
528 \
529 \
530 \
531 \
533 \
534 \
535 \
536 \
537 \
538 \
539 \
540 \
541 \
542 \
543 \
544 \
545 \
546 \
547 \
548 \
549 \
550 \
551 \
552 \
553 \
554 \
555 \
556 \
557 \
558 \
559 \
561 \
562 \
563 \
564 \
565 \
566 \
567 \
568 \
569 \
570 \
571 \
572 \
573 \
574 \
575 \
576 \
577 \
578 \
579 \
580 \
581 \
582 \
583 \
584 \
585 \
586 \
590 \
591 \
592 \
593 \
594 \
595 \
596 \
597 \
598 \
599 \
600 \
601 \
602 \
603 \
604 \
606 \
607 \
608 \
609 \
610 \
611 \
612 \
613 \
614 \
618 \
619 \
620 \
621 \
622 \
623 \
624 \
625 \
626 \
627 \
628 \
629 \
630 \
631 \
632 \
633 \
634 \
635 \
636 \
637 \
638 \
639 \
640 \
641 \
642 \
646 \
647 \
648 \
649 \
650 \
651 \
652 \
653 \
654 \
655 \
656 \
657 \
658 \
659 \
660 \
661 \
662 \
663 \
664 \
665 \
666 \
667 \
668 \
669 \
670 \
675 \
676 \
677 \
678 \
679 \
680 \
681 \
682 \
683 \
684 \
685 \
686 \
687 \
688 \
689 \
690 \
691 \
692 \
693 \
694 \
695 \
696 \
697 \
703 \
704 \
705 \
706 \
707 \
708 \
709 \
710 \
711 \
712 \
713 \
714 \
715 \
716 \
717 \
718 \
719 \
720 \
721 \
722 \
723 \
724 \
725 \
732 \
733 \
734 \
735 \
736 \
737 \
738 \
739 \
740 \
741 \
742 \
743 \
744 \
745 \
746 \
747 \
748 \
749 \
750 \
751 \
752 \
762 \
763 \
764 \
765 \
766 \
767 \
768 \
769 \
770 \
771 \
772 \
773 \
774 \
775 \
776 \
777 

