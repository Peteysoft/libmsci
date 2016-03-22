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

# control as in control experiment:
RAW = test.a
OPTIONS0 = -W 40 -k 400

TRN_CONTR = sq_control$(VER).$(MBTRAIN_SFX)

#we preppend this to the names of all the expendable files involved in 
#multi-borders classification
#to make it easy to clean up:
BASE = test_multi$(VER)/

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

SCRIPT = $(BASE)script$(VER).sh
OUTPUT = $(BASE)borders$(VER)
CLS_CONTR = $(BASE)cls_control$(VER)$(MBCONTROL_SFX)

TRAIN_TIME = $(BASE)train_time$(VER).txt
CLASS_TIME = class_time$(VER).txt

all: $(TEST).m$(VER).cls $(RAW).cls $(TRAIN_TIME) $(CLASS_TIME) raw_time.txt results$(VER).txt
	@echo "Truth vs. AGF (raw):"
	cls_comp_stats $(TEST).cls $(RAW)
	@echo "Time for raw AGF:"
	@more raw_time.txt
	@echo
	@echo "Truth vs. multi-class:"
	cls_comp_stats $(TEST).cls $(TEST).m$(VER)
	@echo "AGF vs. multi-class:"
	cls_comp_stats -b $(RAW).cls $(TEST).m$(VER)
	@echo "AGF vs. multi-class conditional probabilities:"
	agf_correlate $(RAW).con $(TEST).m$(VER).con
	@echo
	@echo "Training time:"
	@more $(TRAIN_TIME)
	@echo
	@echo "Classification time:"
	@more $(CLASS_TIME)

clean:
	rm -f $(BASE)*
	rm -f test.m$(VER).*
	rm -f results$(VER).txt
	rm -f $(CLASS_TIME)
	rm -f classerr$(VER).txt

clean_database:
	rm -f $(DATABASE).cls
	rm -f $(TEST).cls $(TEST).vec
	rm -f $(TRAIN).cls $(TRAIN).vec $(TRAIN).std
	rm -f $(RAW).*
	rm -f results0.txt
	rm -f raw_time.txt
	rm -f rawerr.txt

$(TEST).m$(VER).cls $(TEST).m$(VER).con results$(VER).txt $(CLASS_TIME): $(CLS_CONTR) $(TEST).vec
	( time ( classify_m $(CLS_CONTR) $(TEST).vec $(TEST).m$(VER) > results$(VER).txt 2> classerr$(VER).txt ) ) 2> $(CLASS_TIME)

$(CLS_CONTR) $(TRAIN_TIME): $(TRAIN).cls $(TRAIN).vec $(TRN_CONTR)
	if [ ! -d $(BASE) ]; then mkdir $(BASE); fi
	multi_borders $(OPTIONS$(VER)) $(TRN_CONTR) $(TRAIN) $(OUTPUT) \
		$(CLS_CONTR) > $(SCRIPT)
	chmod u+x $(SCRIPT)
	( time ( ./$(SCRIPT) 2> $(BASE)/trainerr.txt ) ) 2> $(TRAIN_TIME)

$(RAW).cls $(RAW).vec results0.txt raw_time.txt: $(TRAIN).vec $(TRAIN).cls $(TEST).vec
	( time ( agf $(OPTIONS0) classify $(TRAIN) $(TEST).vec $(RAW) > results0.txt 2> rawerr.txt ) ) 2> raw_time.txt

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
LSVMTRAIN_SFX = .lst		#lIBsVM tRAINING CONTROL FILE
LSVMCONTROL_SFX = .lsc		#lIBsVM CLASSIFICATION cONTROL FILE 

#file names:

#base name for training data
SVMTRAIN = svmtrain
#base name for test data
SVMTEST = svmtest

#prepend to SVM stuff:
SVMBASE = sq_svm$(VER)/

#output control file
SVMCONTROL = $(SVMBASE)control$(VER)$(LSVMCONTROL_SFX)
#training script
SVMSCRIPT = $(SVMBASE)svmscript$(VER).sh
#base name for model data
SVMMODEL = $(SVMBASE)svmmodel$(VER)
#output classifications
SVMOUT = svmtest.svm$(VER)

#output from time command:
SVMTRAINTIME=$(SVMBASE)svmtraintime$(VER).txt
SVMCLASSTIME=svmclasstime$(VER).txt

# "native" SVM (one-against-one, isn't passed through AGF multi-borders):
SVMRAW = svmtest0
SVMRAWMODEL = svmmodel0.svm
RAWSVMTRAINTIME=rawsvmtraintime.txt
RAWSVMCLASSTIME=rawsvmclasstime.txt

#"accelerated" version:
ACCBASE = sq_acc$(VER)/
ACCSCRIPT = $(ACCBASE)svm_acc$(VER).sh
ACCCONTROL = $(ACCBASE)svm_acc$(VER).txt
ACCMODEL = $(ACCBASE)svmaccel$(VER)
ACCTRAINTIME=$(ACCBASE)acctraintime$(VER).txt
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
	rm -f $(SVMRAWMODEL)
	rm -f $(SVMTEST).vec $(SVMTEST).cls $(SVMTEST).svm
	rm -f $(SVMRAW).cls $(SVMRAW).con $(SVMRAW).svmout
	rm -f $(RAWSVMTRAINTIME) $(RAWSVMCLASSTIME)
	rm -f $(SVMBASE)*
	rm -f $(SVMOUT).cls $(SVMOUT).con $(SVMOUT).svmout
	rm -f results.svm$(VER).txt
	rm -f $(SVMCLASSTIME)
	rm -f svmerr$(VER).txt
	rm -f $(ACCBASE)*
	rm -f $(ACCOUT).cls $(ACCOUT).con $(ACCOUT).svmout
	rm -f results.acc$(VER).txt
	rm -f $(ACCCLASSTIME)
	rm -f accerr$(VER).txt

#perform classifications with classify_m calling svm-predict:
$(SVMOUT).cls $(SVMOUT).con $(SVMCLASSTIME): $(SVMCONTROL)
	( time ( classify_m -M -O "svm-predict -b 1" $(SVMCONTROL) $(SVMTEST).svm \
		$(SVMOUT).svmout > results.svm$(VER).txt 2> svmerr$(VER).txt ) ) 2> $(SVMCLASSTIME)
	svmout2agf $(SVMOUT).svmout $(SVMOUT)

#create script, model control file and perform training with svm-train:
$(SVMCONTROL) $(SVMSCRIPT) $(SVMTRAINTIME): $(SVMTRAIN).svm $(SVMTEST).svm $(TRN_CONTR)
	if [ ! -d $(SVMBASE) ]; then mkdir $(SVMBASE); fi
	multi_borders -M -- svm-train -+ "-b 1 $(SVMTRAINOPT)" $(TRN_CONTR) \
		$(SVMTRAIN).svm $(SVMMODEL) $(SVMCONTROL) > $(SVMSCRIPT)
	chmod u+x $(SVMSCRIPT)
	( time ( ./$(SVMSCRIPT) 2> $(SVMBASE)svmtrainerr.txt ) ) 2> $(SVMTRAINTIME)

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

#perform classifications with classify_m:
$(ACCOUT).cls $(ACCOUT).con $(ACCCLASSTIME): $(ACCCONTROL) $(SVMTEST).vec
	( time ( classify_m $(ACCCONTROL) $(SVMTEST).vec $(ACCOUT) \
		> results.acc$(VER).txt 2> accerr$(VER).txt ) ) 2> $(ACCCLASSTIME)

#create script, model control file and perform training with class_borders 
#calling svm-predict:
$(ACCSCRIPT) $(ACCCONTROL) $(ACCTRAINTIME): $(SVMCONTROL) $(SVMTRAIN).svm
	if [ ! -d $(ACCBASE) ]; then mkdir $(ACCBASE); fi
	multi_borders -M -O "svm-predict -b 1" $(SVMCONTROL) $(SVMTRAIN).svm \
		$(ACCMODEL) $(ACCCONTROL) > $(ACCSCRIPT)
	chmod u+x $(ACCSCRIPT)
	( time ./$(ACCSCRIPT) ) 2> $(ACCTRAINTIME)

