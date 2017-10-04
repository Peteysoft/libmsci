SHELL=/bin/bash
TIME=/bin/time

DATAPATH=../statlog/data/

VER=shuttle

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

#options for different datasets:
SVMOPT_sat=-c 50 -g 0.1
SVMOPT_segment=-c 100 -g 0.1 
SVMOPT_pendigits=-c 50 -g 0.01
SVMOPT_letter=-c 50 -h 0
SVMOPT_humidity=-h 0 -c 50

OPTIONS=$(SVMOPT_$(VER))

BASE=$(WORKDIR)/$(VER)

TRAIN=$(BASE).trn
TEST=$(BASE).tst

BASEORTHO=$(BASE).ortho
BASE1V1=$(BASE).1v1

#control files:
#training:
CONTROL1V1 = $(WORKDIR)/onevsone.mbt
CONTROLORTHO=$(WORKDIR)/orthogonal.mbt
#classification:
MODELORTHO=$(BASEORTHO).mbc
MODEL1V1=$(BASE1V1).mbc

#timing:
#training time:
#for orthogonal model:
TM_TRAIN_ORTHO=$(WORKDIR)/tm_train_ortho.txt
#for 1 v 1 model:
TM_TRAIN_1V1=$(WORKDIR)/tm_train_1v1.txt
#classification time:
#for orthogonal model:
TM_CLASS_ORTHO=$(WORKDIR)/tm_class_ortho.txt
#for 1 v 1 model:
TM_CLASS_1V1=$(WORKDIR)/tm_class_1v1.txt

#training 1 vs. 1 model:
$(MODEL1V1): $(TRAIN).cls $(TRAIN).vec $(CONTROL1V1)
	$(TIME) -o $(TM_TRAIN_1V1) \
		multi_borders -M -- svm-train -+ "-b 1 $(SVMOPT)" \
		$(CONTROL1V1) $(TRAIN) $(MODEL1V1) $(BASE1V1) $(MODEL1V1)

#training orthogonal model:
$(MODELORTHO): $(TRAIN).cls $(TRAIN).vec $(CONTROLORTHO)
	$(TIME) -o $(TM_TRAIN_ORTHO) \
		multi_borders -M -- svm-train -+ "-b 1 $(SVMOPT)" \
		$(CONTROLORTHO) $(TRAIN) $(CONTROLORTHO) $(BASEORTHO) $(MODELORTHO)

#create control files:
$(CONTROL1V1):
	print_control -Q 6 $(NCLS) > $(CONTROL1V1)

$(CONTROLORTHO):
	print_control -Q 10 $(NCLS) > $(CONTROLORTHO)

#divide into test and training:
$(TRAIN).svm $(TEST).svm: $(BASE).svm
	agf_preprocess -f $(FRAC) $(BASE) $(TRAIN) $(TEST)
