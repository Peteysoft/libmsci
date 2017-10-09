SHELL=/bin/bash
TIME=/bin/time

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
NCLS_segment=7
NCLS_sat=6
NCLS_poker=10
NCLS_mnist=10
NCLS_letter=26
NCLS_pendigits=10
NCLS_usps=10
#this one actually sets the number of classes:
NCLS_humidity=8

NCLS=$(NCLS_$(VER))

#options for AGF-borders:
AGFOPT_shuttle=-W 1.5 -k 200 -s 250
AGFOPT_covtype=
AGFOPT_segment=-W 3 -s 100 -k 150
AGFOPT_sat=-W 5 -s 200 -k 100
AGFOPT_poker=-k 200 -W 2
AGFOPT_mnist=-k 400 -W 25
AGFOPT_vehicle=-W 5 -s 100 -k 100
AGFOPT_letter=-W 3 -k 200 -s 1000
AGFOPT_pendigits=-k 200 -W 2 -s 200
AGFOPT_usps=-k 100 -W 1 -s 250
AGFOPT_ijcnn1=-W 10 -k 400 -s 500
AGFOPT_humidity=-k 400 -W 10 -s 200
AGFOPT_mnist=-k 200 -W 1 -s 500

OPTIONS_agf=$(AGFOPT_$(VER))

#options for different datasets:
SVMOPT_sat=-c 50 -g 0.1
SVMOPT_segment=-c 100 -g 0.1 
SVMOPT_pendigits=-c 50 -g 0.01
SVMOPT_letter=-c 50 -h 0
SVMOPT_humidity=-h 0 -c 50

OPTIONS_svm=$(SVMOPT_$(VER))

OPTIONS=$(OPTIONS_$(METHOD))

WORKDIR=$(VER)_$(METHOD)
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


#command for training:
TRAINCOMMAND_svm=multi_borders -M -- svm-train -+ "-b 1 $(OPTIONS)" \
		$(CONTROL) $(TRAIN).svm $(BASE).$(TYPE) $(MODEL)
TRAINCOMMAND_agf=multi_borders $(OPTIONS) $(CONTROL) $(TRAIN) $(BASE).$(TYPE) $(MODEL)
TRAINCOMMAND=$(TRAINCOMMAND_$(METHOD))

TESTCOMMAND_svm=classify_m -Q $(CODE) -M -O "svm-predict -b 1" \
		$(MODEL) $(TEST).svm $(RESULT).svmout
TESTCOMMAND_agf=classify_m -Q $(CODE) $(MODEL) $(TEST).vec $(RESULT) > $(RESULT).txt
TESTCOMMAND=$(TESTCOMMAND_$(METHOD))

FCONV_svm=svmout2agf $(RESULT).svmout $(RESULT) > $(RESULT).txt
FCONV=$(FCONV_$(METHOD))

#timing:
#training time:
TM_TRAIN=$(BASE).$(TYPE).tm
#classification time:
TM_CLASS=$(BASE).$(TYPE).tm


all:
	make TYPE=1v1 $(BASE).1v1.txt
	make TYPE=ortho $(BASE).ortho.txt
	echo "1 vs. 1 results:"
	more $(BASE).1v1.tm
	cls_comp_stats $(TEST).cls $(BASE).1v1
	echo "orthogonal results:"
	more $(BASE).ortho.tm
	cls_comp_stats $(TEST).cls $(BASE).ortho

clean:
	rm -f $(WORKDIR)/*

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
	agf_preprocess -zf $(FRAC) $(DATAPATH)/$(VER) $(TRAIN).0 $(TEST).0

#make the working directory:
$(WORKDIR)/timestamp:
	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	date > $(WORKDIR)/timestamp

$(DATAPATH)/$(VER).vec:
	make -C $(DATAPATH)
