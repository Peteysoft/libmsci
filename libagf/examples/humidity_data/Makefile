VER = 2

SHELL = /bin/bash

#data set:
DATABASE = est_bt_cld_land_z0.00
ORDFILE = train_SQ36s.dat

#base name of training data:
TRAIN = train
#base name of test data:
TEST = test

#fraction of data to use for testing:
FRAC = 0.1

# control as in control experiment:
CONTROL = test.a
OPTIONS0 = -W 40 -k 400

TRN_CONTR = sq_control$(VER).txt
OPTIONS1 = -s 100 -W 40 -k 450
OPTIONS3 = -s 100 -W 40 -k 400
OPTIONS4 = -s 125 -W 25 -k 300
OPTIONS5 = -s 100 -W 40 -k 400
OPTIONS6 = $(OPTIONS5)
OPTIONS10 = -W 10 -k 400
SCRIPT = script$(VER).sh
OUTPUT = borders$(VER)
CLS_CONTR = cls_control$(VER).txt

TRAIN_TIME = train_time$(VER).txt
CLASS_TIME = class_time$(VER).txt

MIN = 0.00007
MAX = 0.001

NCLS = 8

all: $(TEST).m$(VER).cls $(CONTROL).cls $(TRAIN_TIME) $(CLASS_TIME) raw_time.txt
	@echo "Truth vs. multi-class:"
	cls_comp_stats $(TEST).cls $(TEST).m$(VER)
	@echo "AGF vs. multi-class:"
	cls_comp_stats $(CONTROL).cls $(TEST).m$(VER)
	@echo "AGF vs. multi-class conditional probabilities:"
	agf_correlate $(CONTROL).con $(TEST).m$(VER).con
	@echo
	@echo "Training time:"
	@more $(TRAIN_TIME)
	@echo
	@echo "Classification time:"
	@more $(CLASS_TIME)
	@echo
	@echo "Time for raw AGF:"
	@more raw_time.txt

$(TEST).m$(VER).cls $(CLASS_TIME): $(CLS_CONTR) $(TEST).vec
	( time ( classify_m $(CLS_CONTR) $(TEST).vec $(TEST).m$(VER) > results$(VER).txt 2> classerr.txt ) ) 2> $(CLASS_TIME)

clean:
	rm -f $(SCRIPT)
	rm -f $(OUTPUT)*.b?d
	rm -f $(OUTPUT).std
	rm -f $(CLS_CONTR)
	rm -f $(TEST).m$(VER).cls $(TEST).m$(VER).vec $(TEST).m$(VER).con
	rm -f $(TRAIN_TIME)
	rm -f $(CLASS_TIME)
	rm -f classerr.txt
	rm -f rawerr.txt
	rm -f results$(VER).txt
	rm -f trainerr.txt 

clean_database:
	rm -f $(DATABASE).cls
	rm -f $(TEST).cls $(TEST).vec
	rm -f $(TRAIN).cls $(TRAIN).vec $(TRAIN).std
	rm -f $(CONTROL).*
	rm -f results0.txt
	rm -f raw_time.txt

$(CLS_CONTR) $(TRAIN_TIME): $(TRAIN).cls $(TRAIN).vec $(TRN_CONTR)
	multi_borders $(OPTIONS$(VER)) $(TRN_CONTR) $(TRAIN) $(OUTPUT) \
		$(CLS_CONTR) > $(SCRIPT)
	chmod u+x $(SCRIPT)
	( time ( ./$(SCRIPT) 2> trainerr.txt ) ) 2> $(TRAIN_TIME)

$(TRAIN).vec $(TRAIN).cls $(TEST).vec $(TEST).cls: $(DATABASE).cls $(DATABASE).vec
	agf_preprocess -R -n -f $(FRAC) $(DATABASE) $(TRAIN) $(TEST)

$(DATABASE).cls: $(ORDFILE)
	float_to_class -q $(NCLS) -g $(ORDFILE) $(DATABASE).cls $(MIN) $(MAX)

$(CONTROL).cls $(CONTROL).vec raw_time.txt: $(TRAIN).vec $(TRAIN).cls $(TEST).vec
	( time ( agf $(OPTIONS0) classify $(TRAIN) $(TEST).vec $(CONTROL) > results0.txt 2> rawerr.txt ) ) 2> raw_time.txt

#################### SVM stuff ############################

#fraction of original set to use for this experiment:
SVMFRAC = 0.2

#options to pass to svm-train:
SVMTRAINOPT = -h 0 -c 50

#file names:

#base name for training data
SVMTRAIN = svmtrain
#base name for test data
SVMTEST = svmtest
#output control file
SVMCONTROL = svmcontrol$(VER).txt
#training script
SVMSCRIPT = svmscript$(VER).sh
#base name for model data
SVMMODEL = svmmodel$(VER)
#output classifications
SVMOUT = svmtest.svm$(VER)

#output from time command:
SVMTRAINTIME=svmtraintime$(VER).txt
SVMCLASSTIME=svmclasstime$(VER).txt

#"accelerated" version:
ACCSCRIPT = svm_acc$(VER).sh
ACCCONTROL = svm_acc$(VER).txt
ACCMODEL = svmaccel$(VER)
ACCOUT = svmtest.acc$(VER)
ACCTRAINTIME=acctraintime$(VER).txt
ACCCLASSTIME=accclasstime$(VER).txt

# "native" SVM (one-against-one, isn't passed through AGF multi-borders):
SVMRAW = svmtest0
SVMRAWMODEL = svmmodel0.svm
RAWSVMTRAINTIME=rawsvmtraintime.txt
RAWSVMCLASSTIME=rawsvmclasstime.txt

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
	rm -f $(SVMTRAIN).vec $(SVMTRAIN).cls $(SVMTRAIN).svm
	rm -f $(SVMTEST).vec $(SVMTEST).cls $(SVMTEST).svm
	rm -f $(SVMOUT).cls $(SVMOUT).con $(SVMOUT).svmout
	rm -f $(ACCOUT).cls $(ACCOUT).con $(ACCOUT).svmout
	rm -f $(SVMRAW).cls $(SVMRAW).con $(SVMRAW).svmout
	rm -f $(SVMCONTROL) $(SVMSCRIPT)
	rm -f $(ACCCONTROL) $(ACCSCRIPT)
	rm -f $(SVMMODEL)*.b?d
	rm -f $(SVMRAWMODEL)
	rm -f $(ACCMODEL)*.b?d
	rm -f $(SVMTRAINTIME) $(SVMCLASSTIME)
	rm -f $(ACCTRAINTIME) $(ACCCLASSTIME)
	rm -f $(RAWSVMTRAINTIME) $(RAWSVMCLASSTIME)
	rm -f results.acc$(VER).txt
	rm -f results.svm$(VER).txt
	rm -f svmerr.txt svmtrainerr.txt
	rm -f accerr.txt
	rm -f therest.vec therest.cls

#because svm is effing slow we use only a fraction of the data:
$(SVMTRAIN).svm $(SVMTEST).svm $(SVMTRAIN).vec $(SVMTRAIN).cls $(SVMTEST).vec $(SVMTEST).cls: \
		$(TRAIN).vec $(TRAIN).cls
	agf_preprocess -z -f $(SVMFRAC) $(TRAIN) therest svmall
	agf_preprocess -R -f $(FRAC) svmall $(SVMTRAIN) $(SVMTEST)
	agf2ascii -M $(SVMTRAIN) > $(SVMTRAIN).svm
	agf2ascii -M $(SVMTEST) > $(SVMTEST).svm

#create script, model control file and perform training with svm-train:
$(SVMCONTROL) $(SVMSCRIPT) $(SVMTRAINTIME): $(SVMTRAIN).svm $(SVMTEST).svm $(TRN_CONTR)
	multi_borders -M -- svm-train -+ "-b 1 $(SVMTRAINOPT)" $(TRN_CONTR) \
		$(SVMTRAIN).svm $(SVMMODEL) $(SVMCONTROL) > $(SVMSCRIPT)
	chmod u+x $(SVMSCRIPT)
	( time ( ./$(SVMSCRIPT) 2> svmtrainerr.txt ) ) 2> $(SVMTRAINTIME)

#perform classifications with classify_m calling svm-predict:
$(SVMOUT).cls $(SVMOUT).con $(SVMCLASSTIME): $(SVMCONTROL)
	( time ( classify_m -M -O "svm-predict -b 1" $(SVMCONTROL) $(SVMTEST).svm \
		$(SVMOUT).svmout > results.svm$(VER).txt 2> svmerr.txt ) ) 2> $(SVMCLASSTIME)
	svmout2agf $(SVMOUT).svmout $(SVMOUT)

#"accelerated" SVM:

#create script, model control file and perform training with class_borders 
#calling svm-predict:
$(ACCSCRIPT) $(ACCCONTROL) $(ACCTRAINTIME): $(SVMCONTROL) $(SVMTRAIN).svm
	multi_borders -M -O "svm-predict -b 1" $(SVMCONTROL) $(SVMTRAIN).svm \
		$(ACCMODEL) $(ACCCONTROL) > $(ACCSCRIPT)
	chmod u+x $(ACCSCRIPT)
	( time ./$(ACCSCRIPT) ) 2> $(ACCTRAINTIME)

#perform classifications with classify_m:
$(ACCOUT).cls $(ACCOUT).con $(ACCCLASSTIME): $(ACCCONTROL) $(SVMTEST).vec
	( time ( classify_m $(ACCCONTROL) $(SVMTEST).vec $(ACCOUT) \
		> results.acc$(VER).txt 2> accerr.txt ) ) 2> $(ACCCLASSTIME)

# native SVM multi-class:
#svm training with svm-train:
$(SVMRAWMODEL) $(RAWSVMTRAINTIME): $(SVMTRAIN).svm
	( time svm-train -b 1 $(SVMTRAINOPT) $(SVMTRAIN).svm $(SVMRAWMODEL) ) \
		2> $(RAWSVMTRAINTIME)

#classifications with svm-predict:
$(SVMRAW).cls $(SVMRAW).con $(RAWSVMCLASSTIME): $(SVMRAWMODEL) $(SVMTEST).svm
	(time svm-predict -b 1 $(SVMTEST).svm $(SVMRAWMODEL) $(SVMRAW).svmout ) 2> $(RAWSVMCLASSTIME)
	svmout2agf $(SVMRAW).svmout $(SVMRAW)

