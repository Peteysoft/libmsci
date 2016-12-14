VER = 2

SHELL = /bin/bash

#for dividing up the continuum data:
MIN = 0.00007
MAX = 0.001

#number of classes:
NCLS = 8

#fraction of data to use for testing:
FRAC = 0.1

#data set:
DATABASE = est_bt_cld_land_z0.00
ORDFILE = train_SQ36s.dat

#base name of training data:
TRAIN = train
#base name of test data:
TEST = test

#file extensions:
MBTRAIN_SFX = .mbt		#mULTI-bORDERS tRAINING CONTROL FILE
MBCONTROL_SFX = .mbc		#mULTI-bORDERS CLASSIFICATION cONTROL FILE
TRN_CONTR = sq_control$(VER)$(MBTRAIN_SFX)

# control as in control experiment:
RAWBASE = agfraw
RAWRESULT = $(RAWBASE)
RAWOPT = -W 40 -k 400
RAWTIME = $(RAWBASE).tm

#we preppend this to the names of all the expendable files involved in 
#multi-borders classification
#to make it easy to clean up:
BASE = agf_mb$(VER)

#really need default for this...
OPTIONS1 = -s 100 -W 40 -k 450
OPTIONS2a = -s 100 -W 25 -k 400
OPTIONS2b = -s 100 -W 25 -k 400
OPTIONS3 = -s 100 -W 40 -k 400
OPTIONS4 = -s 125 -W 25 -k 300
OPTIONS5 = -s 100 -W 40 -k 400
OPTIONS6 = $(OPTIONS5)
OPTIONS7 = -s 100 -W 40 -k 450
OPTIONS7a = -s 100 -W 40 -k 450
OPTIONS8 = -s 100 -W 40 -k 450
OPTIONS9 = -s 100 -W 5 -k 250
OPTIONS10 = -W 10 -k 400
OPTIONS11 = -s 100 -W 20 -k 450
OPTIONS11a = -s 100 -W 25 -k 400
OPTIONS12 = -s 100 -W 25 -k 400
OPTIONS12a = -s 100 -W 25 -k 400
OPTIONS12b = -s 200 -W 25 -k 400

#directory containing the AGF multi-borders model:
MODELDIR = $(BASE)
#base names for the binary classification models:
MODELBASE =$(MODELDIR)/$(BASE)

#control file for making classifications:
CLS_CONTR = $(MODELDIR)/$(BASE)$(MBCONTROL_SFX)

#time it takes to train model
TRAIN_TIME = $(MODELDIR)/$(BASE).tm

all: $(BASE).cls $(RAWBASE).cls $(TRAIN_TIME) $(BASE).tm $(RAWBASE).tm $(BASE).txt
	@echo "Truth vs. AGF (raw):"
	cls_comp_stats $(TEST).cls $(RAWBASE)
	@echo "Time for raw AGF:"
	@more $(RAWBASE).tm
	@echo
	@echo "Truth vs. multi-class:"
	cls_comp_stats $(TEST).cls $(BASE)
	@echo "AGF vs. multi-class:"
	cls_comp_stats -b $(RAWBASE).cls $(BASE)
	@echo "AGF vs. multi-class conditional probabilities:"
	agf_correlate $(RAWBASE).con $(BASE).con
	@echo
	@echo "Training time:"
	@more $(TRAIN_TIME)
	@echo
	@echo "Classification time:"
	@more $(BASE).tm

clean_all: clean clean_database clean_raw

clean:
	rm -f $(MODELDIR)/*
	if [ -d $(MODELDIR) ]; then rmdir $(MODELDIR); fi
	rm -f $(BASE)/*

clean_database:
	rm -f $(DATABASE).cls
	rm -f $(TEST).cls $(TEST).vec
	rm -f $(TRAIN).cls $(TRAIN).vec $(TRAIN).std

clean_raw:
	rm -f $(RAWBASE).*

$(BASE).cls $(BASE).con $(BASE).txt $(BASE).tm: $(CLS_CONTR) $(TEST).vec
	( time ( classify_m $(CLS_CONTR) $(TEST).vec $(BASE) > $(BASE).txt 2> $(BASE).log ) ) 2> $(BASE).tm

$(CLS_CONTR) $(TRAIN_TIME): $(TRAIN).cls $(TRAIN).vec $(TRN_CONTR)
	if [ ! -d $(MODELDIR)/ ]; then mkdir $(MODELDIR)/; fi
	( time (multi_borders $(OPTIONS$(VER)) $(TRN_CONTR) $(TRAIN) $(MODELBASE) \
		$(CLS_CONTR) 2> $(MODELBASE).log ) ) 2> $(TRAIN_TIME)

$(RAWBASE).cls $(RAWBASE).vec $(RAWBASE).txt $(RAWBASE).tm: $(TRAIN).vec $(TRAIN).cls $(TEST).vec
	( time ( agf $(RAWOPT) classify $(TRAIN) $(TEST).vec $(RAWBASE) > $(RAWBASE).txt 2> $(RAWBASE).log ) ) 2> $(RAWBASE).tm

$(TRAIN).vec $(TRAIN).cls $(TEST).vec $(TEST).cls: $(DATABASE).cls $(DATABASE).vec
	agf_preprocess -R -n -f $(FRAC) $(DATABASE) $(TRAIN) $(TEST)

$(DATABASE).cls: $(ORDFILE)
	float_to_class -q $(NCLS) -g $(ORDFILE) $(DATABASE).cls $(MIN) $(MAX)

#################### SVM stuff ############################

#fraction of original set to use for this experiment:
SVMFRAC = 0.2

#options to pass to svm-train:
SVMTRAINOPT = -h 0 -c 50

#file extensions:
LSVMTRAIN_SFX = .svt		#svM tRAINING CONTROL FILE
LSVMCONTROL_SFX = .svc		#svM CLASSIFICATION cONTROL FILE 

#file names:

#base name for training data
SVMTRAIN = svmtrain
#base name for test data
SVMTEST = svmtest

#prepend to SVM stuff:
SVMBASE = svm_multi$(VER)
SVMDIR = $(SVMBASE)

#output control file
SVMCONTROL = $(SVMDIR)/$(SVMBASE)(LSVMCONTROL_SFX)
#base name for model data
SVMMODEL = $(SVMDIR)/$(SVMBASE)
#output classifications
SVMOUT = $(SVMBASE)

#output from time command:
SVMTRAINTIME=$(SVMDIR)/$(SVMBASE).tm
SVMCLASSTIME=$(SVMBASE).tm

# "native" SVM (one-against-one, isn't passed through AGF multi-borders):
SVMRAW = svm
SVMRAWMODEL = $(SVMRAW).svmmod
RAWSVMTRAINTIME=$(SVMRAW).trn.tm
RAWSVMCLASSTIME=$(SVMRAW).tst.tm

#"accelerated" version:
ACCBASE = acc_mb$(VER)/
ACCDIR = $(ACCBASE)/
ACCSCRIPT = $(ACCDIR)svm_acc$(VER).sh
ACCCONTROL = $(ACCDIR)svm_acc$(VER).txt
ACCMODEL = $(ACCDIR)svmaccel$(VER)
ACCTRAINTIME=$(ACCDIR)acctraintime$(VER).txt
ACCOUT = svmtest.acc$(VER)
ACCCLASSTIME=accclasstime$(VER).txt

svm: $(SVMOUT).cls $(SVMRAW).cls $(ACCOUT).cls \
		$(SVMOUT).con $(SVMRAW).con $(ACCOUT).con \
		$(SVMTRAINTIME) $(SVMCLASSTIME) \
		$(RAWSVMTRAINTIME) $(RAWCLASSTIME) \
		$(ACCTRAINTIME) $(ACCCLASSTIME) 
	@echo "Truth vs. native SVM:"
	cls_comp_stats $(SVMTEST).cls $(SVMRAW)
	@echo "Time for native SVM:"
	@echo "Training:"
	@more $(RAWSVMTRAINTIME)
	@echo
	@echo "Classification:"
	@more $(RAWSVMCLASSTIME)
	@echo
	@echo "Truth vs. SVM multi-borders:"
	cls_comp_stats $(SVMTEST).cls $(SVMOUT)
	@echo "Native SVM vs. SVM multi-borders:"
	cls_comp_stats $(SVMRAW).cls $(SVMOUT)
	@echo "Native SVM vs. SVM multi-borders conditional probabilities:"
	agf_correlate $(SVMRAW).con $(SVMOUT).con
	@echo
	@echo "Training time:"
	@more $(SVMTRAINTIME)
	@echo
	@echo "Classification time:"
	@more $(SVMCLASSTIME)
	@echo
	@echo
	@echo "Truth vs. accelerated SVM:"
	cls_comp_stats $(SVMTEST).cls $(ACCOUT)
	@echo "Native SVM vs. accelerated SVM:"
	cls_comp_stats $(SVMRAW).cls $(ACCOUT)
	@echo "Native SVM vs. accelerated SVM conditional probabilities:"
	agf_correlate $(SVMRAW).con $(ACCOUT).con
	@echo
	@echo "Training time:"
	@more $(ACCTRAINTIME)
	@echo
	@echo "Classification time:"
	@more $(ACCCLASSTIME)

#note to self: to reduce number of files to clean up, should probably
#use some kind of distinctive base name for all files to use as a "template"
svmclean:
	rm -f svmall.cls svmall.vec
	rm -f therest.vec therest.cls
	rm -f $(SVMTRAIN).vec $(SVMTRAIN).cls $(SVMTRAIN).svm
	rm -f $(SVMTEST).vec $(SVMTEST).cls $(SVMTEST).svm
	rm -f $(SVMRAW).*
	rm -f $(SVMBASE).*
	rm -f $(SVMDIR)/*
	if [ -d $(SVMDIR) ]; then rmdir $(SVMDIR); fi

#perform classifications with classify_m calling svm-predict:
$(SVMOUT).cls $(SVMOUT).con $(SVMCLASSTIME): $(SVMCONTROL)
	( time ( classify_m -M -O "svm-predict -b 1" $(SVMCONTROL) $(SVMTEST).svm \
		$(SVMOUT).svmout > $(SVMOUT).txt 2> $(SVMBASE).log ) ) 2> $(SVMCLASSTIME)
	svmout2agf $(SVMOUT).svmout $(SVMOUT)

#create script, model control file and perform training with svm-train:
$(SVMCONTROL) $(SVMTRAINTIME): $(SVMTRAIN).svm $(SVMTEST).svm $(TRN_CONTR)
	if [ ! -d $(SVMDIR) ]; then mkdir $(SVMDIR); fi
	( time ( multi_borders -M -- svm-train -+ "-b 1 $(SVMTRAINOPT)" $(TRN_CONTR) \
		$(SVMTRAIN).svm $(SVMMODEL) $(SVMCONTROL) 2> $(SVMDIR)svmtrainerr.txt ) ) 2> $(SVMTRAINTIME)

# native SVM multi-class:
#classifications with svm-predict:
$(SVMRAW).cls $(SVMRAW).con $(RAWSVMCLASSTIME): $(SVMRAWMODEL) $(SVMTEST).svm
	(time svm-predict -b 1 $(SVMTEST).svm $(SVMRAWMODEL) $(SVMRAW).svmout ) 2> $(RAWSVMCLASSTIME)
	svmout2agf $(SVMRAW).svmout $(SVMRAW)

#svm training with svm-train:
$(SVMRAWMODEL) $(RAWSVMTRAINTIME): $(SVMTRAIN).svm
	( time svm-train -b 1 $(SVMTRAINOPT) $(SVMTRAIN).svm $(SVMRAWMODEL) ) \
		2> $(RAWSVMTRAINTIME)

#because svm is effing slow we use only a fraction of the data:
$(SVMTRAIN).svm $(SVMTEST).svm $(SVMTRAIN).vec $(SVMTRAIN).cls $(SVMTEST).vec $(SVMTEST).cls: \
		$(TRAIN).vec $(TRAIN).cls
	agf_preprocess -z -f $(SVMFRAC) $(TRAIN) therest svmall
	agf_preprocess -R -f $(FRAC) svmall $(SVMTRAIN) $(SVMTEST)
	agf2ascii -M $(SVMTRAIN) > $(SVMTRAIN).svm
	agf2ascii -M $(SVMTEST) > $(SVMTEST).svm

#"accelerated" SVM:
#"accelerated" version:
ACCBASE = acc_mb$(VER)
ACCDIR = $(ACCBASE)
ACCCONTROL = $(ACCDIR)/$(ACCBASE).sbc
ACCMODEL = $(ACCDIR)/$(ACCBASE)
ACCTRAINTIME=$(ACCDIR)/$(ACCBASE).tm
ACCOUT = $(ACCBASE)
ACCCLASSTIME=$(ACCBASE).tm

acc: $(ACCOUT).cls $(ACCOUT).con $(ACCTRAINTIME) $(ACCCLASSTIME) \
	$(SVMTEST).cls $(SVMRAW).cls $(SVMRAW).con
	@echo "Truth vs. accelerated SVM:"
	cls_comp_stats $(SVMTEST).cls $(ACCOUT)
	@echo "Native SVM vs. accelerated SVM:"
	cls_comp_stats $(SVMRAW).cls $(ACCOUT)
	@echo "Native SVM vs. accelerated SVM conditional probabilities:"
	agf_correlate $(SVMRAW).con $(ACCOUT).con
	@echo
	@echo "Training time:"
	@more $(ACCTRAINTIME)
	@echo
	@echo "Classification time:"
	@more $(ACCCLASSTIME)

accclean:
	rm -f $(ACCBASE).*
	rm -f $(ACCDIR)/*
	if [ -d $(ACCDIR) ]; then rmdir $(ACCDIR); fi

#perform classifications with classify_m:
$(ACCOUT).cls $(ACCOUT).con $(ACCCLASSTIME): $(ACCCONTROL) $(SVMTEST).vec
	( time ( classify_m $(ACCCONTROL) $(SVMTEST).vec $(ACCOUT) \
		> $(ACCBASE).txt 2> $(ACCBASE).log ) ) 2> $(ACCCLASSTIME)

#create script, model control file and perform training with class_borders 
#calling svm-predict:
$(ACCSCRIPT) $(ACCCONTROL) $(ACCTRAINTIME): $(SVMCONTROL) $(SVMTRAIN).svm
	if [ ! -d $(ACCDIR) ]; then mkdir $(ACCDIR); fi
	multi_borders -M -O "svm-predict -b 1" $(SVMCONTROL) $(SVMTRAIN).svm \
		$(ACCMODEL) $(ACCCONTROL) > $(ACCSCRIPT)
	chmod u+x $(ACCSCRIPT)
	( time ./$(ACCSCRIPT) ) 2> $(ACCTRAINTIME)

