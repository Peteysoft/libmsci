#minimum number of border samples:
MINSAMPLE=5
#maximum number of border samples:
MAXSAMPLE=200
NTEST=20
NTRIAL=5
#fraction of data used for testing:
FRAC=0.3

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
NCLS_poker=10
NCLS_mnist=10

NCLS=$(NCLS_$(VER))

#special pre-processing:
PRE_heart=-n
PRE_shuttle=-n
PRE_covtype=-n
PRE_segment=-n
PRE_sat=-nS 9
PRE_poker=-n
#PRE_mnist=-nS 20

PRE=$(PRE_$(VER))

#options for AGF-borders:
AGFOPT_heart=-W 5 -k 100
AGFOPT_shuttle=-W 2
AGFOPT_covtype=
AGFOPT_segment=-W 5
AGFOPT_sat=-W 10
AGFOPT_poker=-k 200 -W 2
AGFOPT_mnist=-k 400 -W 25

AGFOPT=$(AGFOPT_$(VER))

#options for svm-train:
SVMOPT_sat=-c 25

SVMOPT=$(SVMOPT_$(VER))

#number of nearest neighbours for knn:
K_heart=-k 11
K_shuttle=-k 1
K_covtype=-k 11
K_segment=-k 3
K_sat=-k 15

K=$(K_$(VER))

#number of border samples (or other options) for "accelerated" SVM model:
NSAMP_segment=-s 20

NSAMP=$(NSAMP_$(VER))

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

MULTICONTROL=onevsone.mbc

SVMMODEL=$(BASE).svmmod
AGFMODEL=$(BASE).agf
ACCELERATOR=$(BASE).accel

NBORDER_AGF=$(WORKDIR)/nborder_agf.txt
NBORDER_SVM=$(WORKDIR)/nborder_svm.txt
NSUPPORT=$(WORKDIR)/nsupport.txt

#need this to normalize test data same as training data, if applicable:
ifneq ($(strip $(PRE)),)
  NORMCLAUSE=-a $(TRAIN).std
endif

#output results:
all: $(TEST).knn.cls $(TEST).agf.cls $(TEST).svm.cls $(TEST).acc.cls $(NBORDER_SVM)
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

clean:
	rm -r $(WORKDIR)/*

#testing:
$(TEST).knn.cls $(TEST).knn.cor: $(TRAIN).vec $(TRAIN).cls $(TEST).vec
	time -o $(WORKDIR)/knn_time.txt knn $K classify $(TRAIN) $(TEST).vec $(TEST).knn

$(TEST).svm.cls $(TEST).svm.cor: $(SVMMODEL) $(TEST).svm
	time -o $(WORKDIR)/svm_time.txt svm-predict -b 1 $(TEST).svm $(SVMMODEL) $(TEST).svm.txt
	svmout2agf $(TEST).svm.txt > $(TEST).svm

$(TEST).agf.cls $(TEST).agf.cor: $(AGFMODEL).mbc $(TEST).vec
	time -o $(WORKDIR)/agf_time.txt classify_m $(AGFMODEL).mbc $(TEST).vec $(TEST).agf > $(TEST).agf.txt

$(TEST).acc.cls $(TEST).acc.cor: $(ACCELERATOR) $(TEST).vec
	time -o $(WORKDIR)/acc_time.txt classify_s $(ACCELERATOR) $(TEST).vec $(TEST).acc > $(TEST).acc.txt

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
$(ACCELERATOR): $(SVMMODEL) $(TRAIN).vec $(TRAIN).cls
	time -o $(WORKDIR)/accel_ttime.txt svm_accelerate $(NSAMP) $(SVMMODEL) $(TRAIN) $(ACCELERATOR)

#svm model:
$(SVMMODEL): $(TRAIN).svm
	time -o $(WORKDIR)/svm_ttime.txt svm-train -b 1 $(SVMOPT) $(TRAIN).svm $(SVMMODEL)

#agf borders model:
$(AGFMODEL).mbc: $(TRAIN).svm
	time -o $(WORKDIR)/agf_ttime.txt print_control -Q 6 $(NCLS) | multi_borders $(AGFOPT) $(TRAIN) $(AGFMODEL) $(AGFMODEL).mbc

#convert to LIBSVM format:
$(TRAIN).svm: $(TRAIN).cls $(TRAIN).vec
	agf2ascii -M $(TRAIN) $(TRAIN).svm

$(TEST).svm: $(TEST).cls $(TEST).vec
	agf2ascii -M $(TEST) $(TEST).svm

#since multiple targets aren't possible with suffix rules (I don't think):
$(TRAIN).vec $(TEST).cls $(TEST).vec: $(TRAIN).cls

#file conversion:
#steps:
#1. convert to libAGF format and copy to working directory
#2. split into test and training, if applicable
#3. pre-processing: 
#     - normalization (optional)
#     - SVD (optional)
#     - make class labels go from 0 to Nc-1
$(WORKDIR)/heart.trn.cls: $(WORKDIR)/%.trn.cls: $(DATADIR)/%.dat.txt
	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	lvq2agf -H $< $(WORKDIR)/$(VER)
	agf_preprocess $(PRE) -Uf $(FRAC) $(WORKDIR)/$(VER) $(TRAIN) $(TEST)

#3rd feature must be removed before normalization 
#since it only takes on one value:
$(WORKDIR)/segment.trn.cls: $(DATADIR)/segment.dat.txt
	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	lvq2agf -H $< $(WORKDIR)/$(VER)
	agf_preprocess $(PRE) -FUf $(FRAC) $(WORKDIR)/$(VER) $(TRAIN) $(TEST) 0 1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18

$(WORKDIR)/sat.trn.cls $(WORKDIR)/shuttle.trn.cls: $(WORKDIR)/%.cls: $(DATADIR)/%.txt
	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	lvq2agf -H $(DATADIR)/$(VER).trn.txt $(TRAIN)
	lvq2agf -H $(DATADIR)/$(VER).tst.txt $(TEST)
	agf_preprocess -U $(PRE) $(TRAIN) $(TRAIN)
	agf_preprocess $(NORMCLAUSE) -U $(TEST) $(TEST)

VEHICLE_FILES=xaa.dat xab.dat xac.dat xad.dat xae.dat xaf.dat xag.dat xah.dat xai.dat
VEHICLE_TEST=xaa.dat
VEHICLE_TRAIN=$(VEHICLE_FILES, $(VEHICLE_TEST)=)

$(WORKDIR)/vehicle.trn.cls: $(addprefix, $(DATADIR)/, $(VEHICLE_FILES))
	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	cat $(addprefix, $(DATADIR)/, $(VEHICLE_TRAIN)) > $(WORKDIR)/vehicle/vehicle.trn.txt
	cp $(addprefix, $(DATADIR)/, VEHICLE_TEST) $(WORKDIR)/vehicle/vehicle.tst.txt
	sed 's/opel/0; s/saab/1; s/bus/2; s/van/3' $(WORKDIR)/$(VER).trn.txt
	sed 's/opel/0; s/saab/1; s/bus/2; s/van/3' $(WORKDIR)//$(VER).tst.txt
	lvq2agf -H $< $(WORKDIR)/$(VER).trn.txt $(TRAIN)
	lvq2agf -H $< $(WORKDIR)/$(VER).tst.txt $(TEST)
	agf_preprocess -U $(PRE) $(TRAIN) $(TRAIN)
	agf_preprocess $(NORMCLAUSE) -U $(TEST) $(TEST)

$(WORKDIR)/covtype.trn.cls: $(WORKDIR)/%.trn.cls: $(DATADIR)/%
	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	svm2agf -E 0 $(DATADIR)/$(VER) $(WORKDIR)/$(VER)
	agf_preprocess $(PRE) -Uf $(FRAC) $(WORKDIR)/$(VER) $(TRAIN) $(TEST)

$(WORKDIR)/mnist.trn.cls $(WORKDIR)/poker.trn.cls: $(WORKDIR)/%.trn.cls: $(DATADIR)/%
	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	svm2agf -E 0 $(DATADIR)/$(VER) $(TRAIN)
	svm2agf -E 0 $(DATADIR)/$(VER).t $(TEST)
	agf_preprocess $(PRE) -U $(TRAIN) $(TRAIN)
	agf_preprocess $(NORMCLAUSE) $(TEST) $(TEST)

