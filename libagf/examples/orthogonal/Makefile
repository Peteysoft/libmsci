SHELL=/bin/bash
#TIME=/bin/time
TIME=/usr/bin/time

DATAPATH=../datasets

VER=shuttle

#which method are we using?
METHOD=svm
#METHOD=agf

#fraction of data to use for training:
FRAC=0.3
#number of trials:
NTRIAL=20

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

OPTIONS=$(OPTIONS_$(METHOD))

WORKDIR=/home/petey/work/orthogonal/$(VER)_$(METHOD)
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

#control file option:
CODE1v1=5
CODEortho=8
CODE=$(CODE$(TYPE))

#preprocessing options:
PRE=-n

#need this to normalize test data same as training data, if applicable:
ifneq ($(strip $(PRE)),)
  NORMFILE=$(TRAIN).std
  NORMCLAUSE=-a $(NORMFILE)
endif

#subsampling option:
CFLAG_shuttle=correction
#keep class sizes constant:
SAMPLING_METHOD_= -C
#correct for smaller classes:
SAMPLING_METHOD_correction=
#there must be a way of doing it without using three levels of indirection
#(and without resorting to 'if'):
SAMPLING_METHOD=$(SAMPLING_METHOD_$(CFLAG_$(VER)))

#command for subsampling:
#takes into account relative class sizes:
SUBSAMPLE=subsample_special
#doesn't:
#SUBSAMPLE=agf_preprocess


#command for training:
TRAINCOMMAND_svm=multi_borders -M -- svm-train -+ "-b 1 $(OPTIONS)" \
		$(CONTROL) $(TRAIN).svm $(BASE).$(TYPE) $(MODEL).0; \
		mbh2mbm -Z $(MODEL).0 $(MODEL)
TRAINCOMMAND_agf=multi_borders $(OPTIONS) $(CONTROL) $(TRAIN) $(BASE).$(TYPE) $(MODEL)
TRAINCOMMAND=$(TRAINCOMMAND_$(METHOD))

#TESTCOMMAND_svm=classify_m -Q $(CODE) -M -O "svm-predict -b 1" \
#		$(MODEL) $(TEST).svm $(RESULT).svmout
TESTCOMMAND_svm=classify_m -Z -Q $(CODE) $(MODEL) $(TEST).vec $(RESULT) > $(RESULT).txt
TESTCOMMAND_agf=classify_m -Q $(CODE) $(MODEL) $(TEST).vec $(RESULT) > $(RESULT).txt
TESTCOMMAND=$(TESTCOMMAND_$(METHOD))

#FCONV_svm=svmout2agf $(RESULT).svmout $(RESULT) > $(RESULT).txt
FCONV=$(FCONV_$(METHOD))

#timing:
#training time:
TM_TRAIN=$(BASE).trn.$(TYPE).tm
#classification time:
TM_CLASS=$(BASE).tst.$(TYPE).tm

#collect statistics:
STATFILE = $(VER).txt


all: $(BASE).1v1.txt $(BASE).ortho.txt
	echo "1 vs. 1 results:"
	more $(BASE).trn.1v1.tm
	more $(BASE).tst.1v1.tm
	cls_comp_stats $(TEST).cls $(BASE).1v1
	validate_probabilities -H $(TEST).cls $(BASE).1v1.txt 
	echo "orthogonal results:"
	more $(BASE).trn.ortho.tm
	more $(BASE).tst.ortho.tm
	cls_comp_stats $(TEST).cls $(BASE).ortho
	validate_probabilities -H $(TEST).cls $(BASE).ortho.txt

clean:
	rm -f $(WORKDIR)/*

$(BASE)$(TYPE).1v1.txt:
	make TYPE=1v1 $(BASE).1v1.txt

$(BASE)$(TYPE).ortho.txt:
	make TYPE=ortho $(BASE).ortho.txt


#how many lines does results file contain so far?
NLINE!=cat $(STATFILE) | wc -l
stats:
	for ((I=$(NLINE); I<$(NTRIAL); I++)); do \
		make clean; \
		make $(STATFILE); \
	done

$(STATFILE): $(BASE).1v1.txt $(BASE).ortho.txt
	grep -o "[0-9]*\.[0-9]*" $(BASE).tst.1v1.tm | sed '1q;d' | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	grep -o "[0-9]*\.[0-9]*" $(BASE).tst.1v1.tm | sed '2q;d' | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	grep -o "[0-9]*\.[0-9]*" $(BASE).tst.ortho.tm | sed '1q;d' | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	grep -o "[0-9]*\.[0-9]*" $(BASE).tst.ortho.tm | sed '2q;d' | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	cls_comp_stats -Hb $(TEST).cls $(BASE).1v1 | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	cls_comp_stats -Hb $(TEST).cls $(BASE).ortho | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	validate_probabilities -BH $(TEST).cls $(BASE).1v1.txt | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	validate_probabilities -BH $(TEST).cls $(BASE).ortho.txt >> $(STATFILE)

#performing classifications:
$(RESULT).txt: $(MODEL) $(TEST).svm
	$(TIME) -o $(TM_CLASS) $(TESTCOMMAND)
	$(FCONV)

#training the model:
$(MODEL): $(TRAIN).svm $(CONTROL)
	$(TIME) -o $(TM_TRAIN) $(TRAINCOMMAND)

#create control files:
$(CONTROL):
	print_control -Q $(CODE) $(NCLS) > $(CONTROL)

#normalize:
$(TRAIN).svm: $(TRAIN).0.vec
	agf_preprocess $(PRE) $(TRAIN).0 $(TRAIN)
	agf2ascii -M $(TRAIN) > $(TRAIN).svm

$(TEST).svm: $(TEST).0.vec
	agf_preprocess $(NORMCLAUSE) $(TEST).0 $(TEST)
	agf2ascii -M $(TEST) > $(TEST).svm

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

