SHELL=/bin/bash
#TIME=/bin/time
TIME=/usr/bin/time

#fraction of data to use for training:
FRAC=0.3
#number of trials:
NTRIAL=20

VER=shuttle

#which method are we using?
#METHOD=svm
METHOD=lin
#METHOD=agf

DATAPATH=../datasets

#WORKDIR=/home/peter/work/orthogonal/$(VER)_$(METHOD)
WORKDIR=work/$(VER)_$(METHOD)
#WORKDIR=$(VER)_$(METHOD)
BASE=$(WORKDIR)/$(VER)

TRAIN=$(BASE).trn
TEST=$(BASE).tst
RESULT=$(BASE).$(TYPE)

#file extensions for results:
EXT_agf=.vec
EXT_svm=.svmout
EXT=$(EXT_$(METHOD))

#control files:
#training:
CONTROL = $(BASE).$(TYPE).mbt

#classification:
MODEL=$(BASE).$(TYPE).mbc

TRAIN_DEP_svm=$(TRAIN).svm $(CONTROL)
TRAIN_DEP_lin=$(TRAIN).svm $(CONTROL)
TRAIN_DEP_agf=$(TRAIN).svm $(CONTROL)
TRAIN_DEP_acc=$(WORKDIR)/$(VER)_svm.$(TYPE).mbc

TRAIN_DEP=$(TRAIN_DEP_$(METHOD))

#number of classes in each dataset:
NCLS_shuttle=7
NCLS_covtype=7
NCLS_covtype_s=7
NCLS_segment=7
NCLS_sat=6
NCLS_poker=10
NCLS_mnist_s=10
NCLS_letter=26
NCLS_pendigits=10
NCLS_usps=10
#this one actually sets the number of classes:
NCLS_humidity=8
NCLS_vehicle=4
NCLS_urban=9

NCLS=$(NCLS_$(VER))

#options for AGF-borders:
AGFOPT_shuttle=-W 1.5 -k 200 -s 250
AGFOPT_covtype=
AGFOPT_segment=-W 3 -s 100 -k 150
AGFOPT_sat=-W 5 -s 200 -k 100
AGFOPT_poker=-k 200 -W 2
AGFOPT_mnist=-k 400 -W 25
AGFOPT_letter=-W 3 -k 200 -s 1000
AGFOPT_pendigits=-k 200 -W 2 -s 200
AGFOPT_usps=-k 100 -W 1 -s 250
AGFOPT_humidity=-k 400 -W 10 -s 200
AGFOPT_mnist=-k 200 -W 1 -s 500
AGFOPT_vehicle=-W 1 -s 200

OPTIONS_agf=$(AGFOPT_$(VER))

#options for different datasets:
SVMOPT_sat=-c 50 -g 0.1
SVMOPT_segment=-c 100 -g 0.1 
SVMOPT_pendigits=-c 50 -g 0.01
SVMOPT_letter=-c 50 -h 0
SVMOPT_humidity=-h 0 -c 50

OPTIONS_svm=$(SVMOPT_$(VER))

#number of border samples (or other options) for "accelerated" SVM model:
NSAMP_segment_1v1=-s 50
NSAMP_vehicle_1v1=-s 50
NSAMP_sat_1v1=-s 200
NSAMP_letter_1v1=-s 75
NSAMP_pendigits_1v1=-s 50
NSAMP_usps_1v1=-s 50
NSAMP_humidity_1v1=-s 200
NSAMP_mnist_1v1=-s 500

NSAMP_segment=-s 100
NSAMP_shuttle=-s 200
NSAMP_vehicle=-s 100
NSAMP_sat=-s 300
NSAMP_letter=-s 150
NSAMP_pendigits=-s 100
NSAMP_usps=-s 100
NSAMP_humidity=-s 300
NSAMP_mnist=-s 500

ifeq ($(strip $(TYPE)),1v1)
  OPTIONS_acc=$(NSAMP_$(VER)_1v1)
else
  OPTIONS_acc=$(NSAMP_$(VER))
endif


#options for lib-linear:
LINOPT_segment=-c 50
LINOPT_shuttle=-c 200
LINOPT_pendigits=-c 200
LINOPT_sat=-c 200
LINOPT_usps=-c 50
LINOPT_vehicle=-c 50
LINOPT_humidity=-c 20
LINOPT_urban=-c 0.05

OPTIONS_lin=$(LINOPT_$(VER))

OPTIONS=$(OPTIONS_$(METHOD))

#options for print_control for orthogonal "non-strict" coding matrices:
PC_PARAM4=7 4
PC_PARAM6=12 6
PC_PARAM7=15 7
PC_PARAM8=17 8
PC_PARAM9=20 9
PC_PARAM10=23 10

PC_PARAM_orthonon=$(PC_PARAM$(NCLS))
PC_PARAM=$(PC_PARAM_$(TYPE))

#control file option:
TRAINCODE1v1=1
TRAINCODE1vR=0
TRAINCODEortho=4
TRAINCODEecc=16
TRAINCODEorthonon=10
TRAINCODE=$(TRAINCODE$(TYPE))

#classify_m option:
TESTCODE1v1=1
TESTCODE1vR=0
TESTCODEortho=4
TESTCODEecc=0
TESTCODEorthonon=0
TESTCODE=$(TESTCODE$(TYPE))

#preprocessing options:
PRE=-n

#need this to normalize test data same as training data, if applicable:
ifneq ($(strip $(PRE)),)
  NORMFILE=$(TRAIN).std
  NORMCLAUSE=-a $(NORMFILE)
endif

#subsampling option:
#CFLAG_shuttle=correction
#keep class sizes constant:
SAMPLING_METHOD_= -C
#correct for smaller classes:
#SAMPLING_METHOD_correction=
#there must be a way of doing it without using three levels of indirection
#(and without resorting to 'if'):
SAMPLING_METHOD=$(SAMPLING_METHOD_$(CFLAG_$(VER)))
#(no need to do this anyway: why the fuck am I doing this?)

#command for subsampling:
#takes into account relative class sizes:
SUBSAMPLE=subsample_special
#doesn't:
#SUBSAMPLE=agf_preprocess


#command for training:
TRAINCOMMAND_svm=multi_borders -M -- svm-train -+ "-b 1 $(OPTIONS)" \
		$(CONTROL) $(TRAIN).svm $(BASE).$(TYPE) $(MODEL).0; \
		mbh2mbm -Z $(MODEL).0 $(MODEL)
TRAINCOMMAND_acc=multi_borders -Z $(OPTIONS) $(WORKDIR)/$(VER)_svm.$(TYPE).mbc $(TRAIN) $(BASE).$(TYPE) $(MODEL)
TRAINCOMMAND_lin=multi_borders -M -- lin-train -+ "-s 0 -B 1 $(OPTIONS)" \
		$(CONTROL) $(TRAIN).svm $(BASE).$(TYPE) $(MODEL)
TRAINCOMMAND_agf=multi_borders $(OPTIONS) $(CONTROL) $(TRAIN) $(BASE).$(TYPE) $(MODEL)
TRAINCOMMAND=$(TRAINCOMMAND_$(METHOD))

#TESTCOMMAND_svm=classify_m -Q $(CODE) -M -O "svm-predict -b 1" \
#		$(MODEL) $(TEST).svm $(RESULT).svmout
TESTCOMMAND_svm=classify_m -Z -Q $(TESTCODE) $(MODEL) $(TEST).vec $(RESULT) > $(RESULT).txt
TESTCOMMAND_acc=classify_m -Q $(TESTCODE) $(MODEL) $(TEST).vec $(RESULT) > $(RESULT).txt
TESTCOMMAND_lin=classify_m -G -Q $(TESTCODE) $(MODEL) $(TEST).vec $(RESULT) > $(RESULT).txt
TESTCOMMAND_agf=classify_m -Q $(TESTCODE) $(MODEL) $(TEST).vec $(RESULT) > $(RESULT).txt
TESTCOMMAND=$(TESTCOMMAND_$(METHOD))

#FCONV_svm=svmout2agf $(RESULT).svmout $(RESULT) > $(RESULT).txt
FCONV=$(FCONV_$(METHOD))

#timing:
#training time:
TM_TRAIN=$(BASE).trn.$(TYPE).tm
#classification time:
TM_CLASS=$(BASE).tst.$(TYPE).tm

#where to store solution times:
SOL_TIME=sol_time_$(METHOD)_$(VER)_$(TYPE).txt

#collect statistics:
STATFILE = $(VER)_$(METHOD).txt


all: $(BASE).1v1.txt $(BASE).1vR.txt $(BASE).ecc.txt $(BASE).ortho.txt $(BASE).orthonon.txt
	echo "1 vs. 1 results:"
	more $(BASE).trn.1v1.tm
	more $(BASE).tst.1v1.tm
	cls_comp_stats $(TEST).cls $(BASE).1v1
	validate_probabilities -H $(TEST).cls $(BASE).1v1.txt 
	echo "1 vs. R results:"
	more $(BASE).trn.1vR.tm
	more $(BASE).tst.1vR.tm
	cls_comp_stats $(TEST).cls $(BASE).1vR
	validate_probabilities -H $(TEST).cls $(BASE).1vR.txt 
	echo "ECC results:"
	more $(BASE).trn.ecc.tm
	more $(BASE).tst.ecc.tm
	cls_comp_stats $(TEST).cls $(BASE).ecc
	validate_probabilities -H $(TEST).cls $(BASE).ecc.txt 
	echo "orthogonal results:"
	more $(BASE).trn.ortho.tm
	more $(BASE).tst.ortho.tm
	cls_comp_stats $(TEST).cls $(BASE).ortho
	validate_probabilities -H $(TEST).cls $(BASE).ortho.txt
	echo "orthogonal \"non-strict\" results:"
	more $(BASE).trn.orthonon.tm
	more $(BASE).tst.orthonon.tm
	cls_comp_stats $(TEST).cls $(BASE).orthonon
	validate_probabilities -H $(TEST).cls $(BASE).orthonon.txt

clean:
	rm -f $(WORKDIR)/*

#relies on TYPE not being defined (kind of ugly hack...):
$(BASE)$(TYPE).1v1.txt:
	make TYPE=1v1 $(BASE).1v1.txt

$(BASE)$(TYPE).1vR.txt:
	make TYPE=1vR $(BASE).1vR.txt

$(BASE)$(TYPE).ecc.txt:
	make TYPE=ecc $(BASE).ecc.txt

$(BASE)$(TYPE).ortho.txt:
	make TYPE=ortho $(BASE).ortho.txt

$(BASE)$(TYPE).orthonon.txt:
	make TYPE=orthonon $(BASE).orthonon.txt

#how many lines does results file contain so far?
NLINE!=cat $(STATFILE) | wc -l
stats:
	for ((I=$(NLINE); I<$(NTRIAL); I++)); do \
		make METHOD=svm clean; \
		make clean; \
		make $(STATFILE); \
	done

$(STATFILE): $(BASE).orthonon.txt $(BASE).1v1.txt $(BASE).1vR.txt $(BASE).ecc.txt $(BASE).ortho.txt
#$(STATFILE): $(BASE).ecc.txt $(BASE).ortho.txt
	grep -o "[0-9]*\.[0-9]*" $(BASE).tst.1v1.tm | sed '1q;d' | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	grep -o "[0-9]*\.[0-9]*" $(BASE).tst.1v1.tm | sed '2q;d' | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	cls_comp_stats -Hb $(TEST).cls $(BASE).1v1 | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	validate_probabilities -BH $(TEST).cls $(BASE).1v1.txt | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	grep -o "[0-9]*\.[0-9]*" $(BASE).tst.1vR.tm | sed '1q;d' | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	grep -o "[0-9]*\.[0-9]*" $(BASE).tst.1vR.tm | sed '2q;d' | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	cls_comp_stats -Hb $(TEST).cls $(BASE).1vR | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	validate_probabilities -BH $(TEST).cls $(BASE).1vR.txt | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	grep -o "[0-9]*\.[0-9]*" $(BASE).tst.ecc.tm | sed '1q;d' | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	grep -o "[0-9]*\.[0-9]*" $(BASE).tst.ecc.tm | sed '2q;d' | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	cls_comp_stats -Hb $(TEST).cls $(BASE).ecc | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	validate_probabilities -BH $(TEST).cls $(BASE).ecc.txt | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	grep -o "[0-9]*\.[0-9]*" $(BASE).tst.ortho.tm | sed '1q;d' | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	grep -o "[0-9]*\.[0-9]*" $(BASE).tst.ortho.tm | sed '2q;d' | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	cls_comp_stats -Hb $(TEST).cls $(BASE).ortho | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	validate_probabilities -BH $(TEST).cls $(BASE).ortho.txt | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	grep -o "[0-9]*\.[0-9]*" $(BASE).tst.orthonon.tm | sed '1q;d' | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	grep -o "[0-9]*\.[0-9]*" $(BASE).tst.orthonon.tm | sed '2q;d' | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	cls_comp_stats -Hb $(TEST).cls $(BASE).orthonon | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	validate_probabilities -BH $(TEST).cls $(BASE).orthonon.txt >> $(STATFILE)

#performing classifications:
$(RESULT).txt: $(MODEL) $(TEST).svm
	export SOL_TIME=$(SOL_TIME); $(TIME) -o $(TM_CLASS) $(TESTCOMMAND)
	$(FCONV)

#training the model:
$(MODEL): $(TRAIN_DEP)
	$(TIME) -o $(TM_TRAIN) $(TRAINCOMMAND)

$(WORKDIR)/$(VER)_svm.$(TYPE).mbc: $(WORKDIR)/timestamp
	make METHOD=svm work/$(VER)_svm/$(VER).$(TYPE).mbc
	cp work/$(VER)_svm/$(VER).$(TYPE).mbc.0 $(WORKDIR)/$(VER)_svm.$(TYPE).mbc
	cp work/$(VER)_svm/$(VER).trn.vec $(TRAIN).vec
	cp work/$(VER)_svm/$(VER).trn.cls $(TRAIN).cls

#create control files:
$(CONTROL):
	print_control -Q $(TRAINCODE) $(NCLS) $(PC_PARAM) > $(CONTROL)

#normalize:
$(TRAIN).svm: $(TRAIN).0.vec
	agf_preprocess $(PRE) $(TRAIN).0 $(TRAIN)
	agf2ascii -M $(TRAIN) > $(TRAIN).svm

ifeq ($(strip $(METHOD)),acc)
  $(TEST).svm:
	make METHOD=svm work/$(VER)_svm/$(VER).tst.svm
	cp work/$(VER)_svm/$(VER).tst.vec $(TEST).vec
	cp work/$(VER)_svm/$(VER).tst.cls $(TEST).cls
	cp work/$(VER)_svm/$(VER).tst.svm $(TEST).svm
else
  $(TEST).svm: $(TEST).0.vec
	agf_preprocess $(NORMCLAUSE) $(TEST).0 $(TEST)
	agf2ascii -M $(TEST) > $(TEST).svm
endif

#divide into test and training:
$(TRAIN).0.vec $(TEST).0.vec: $(DATAPATH)/$(VER).vec $(WORKDIR)/timestamp
	$(SUBSAMPLE) $(SAMPLING_METHOD) -zf $(FRAC) $(DATAPATH)/$(VER) $(TRAIN).0 $(TEST).0

#make the working directory:
$(WORKDIR)/timestamp:
	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	date > $(WORKDIR)/timestamp

$(DATAPATH)/$(VER).vec:
	make -C $(DATAPATH)

#$(DATAPATH)/mnist_s.vec:
#	make -C $(DATAPATH) VER=mnist
#	$(SUBSAMPLE) -zf 0.1 $(DATAPATH)/mnist dum $(DATAPATH)/mnist_s

