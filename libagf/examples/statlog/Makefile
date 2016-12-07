#minimum number of border samples:
MINSAMPLE=5
#maximum number of border samples:
MAXSAMPLE=200
NTEST=10
NTRIAL=10
#fraction of data used for testing:
FRAC=0.2

#VER=heart
#VER=shuttle
VER=covtype

NCLS_heart=2
NCLS_shuttle=7
NCLS_covtype=7

NCLS=$(NCLS_$(VER))

AGFOPT_heart=-W 5 -k 100
AGFOPT_shuttle=-W 2
AGFOPT_covtype=

AGFOPT=$(AGFOPT_$(VER))

SVMOPT1=
SVMOPT2=
SVMOPT3=
SVMOPT4=
SVMOPT5=

SVMOPT=$(SVMOPT_$(VER))

K_heart=11
K_shuttle=1
K_covtype=11

K=$(K_$(VER))

SCRIPTDIR=../../scripts

DATADIR=data/$(VER)

WORKDIR=work/$(VER)

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
	time -o $(WORKDIR)/knn_time.txt knn -k $K classify $(TRAIN) $(TEST).vec $(TEST).knn

$(TEST).svm.cls $(TEST).svm.cor: $(SVMMODEL) $(TEST).svm
	time -o $(WORKDIR)/svm_time.txt svm-predict -b 1 $(TEST).svm $(SVMMODEL) $(TEST).svm.txt
	svmout2agf $(TEST).svm.txt $(TEST).svm

$(TEST).agf.cls $(TEST).agf.cor: $(AGFMODEL).mbc $(TEST).vec
	time -o $(WORKDIR)/agf_time.txt classify_m $(AGFMODEL).mbc $(TEST).vec $(TEST).agf

$(TEST).acc.cls $(TEST).acc.cor: $(ACCELERATOR) $(TEST).vec
	time -o $(WORKDIR)/acc_time.txt classify_s $(ACCELERATOR) $(TEST).vec $(TEST).acc

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

$(NSUPPORT): $(TRAIN).svm
	$(SCRIPTDIR)/svm_nsv.sh -G -q $(NTEST) $(SVMOPT) $(TRAIN) $(NSUPPORT)

#svm-borders "accelerated" svm model
$(ACCELERATOR): $(SVMMODEL) $(TRAIN).vec $(TRAIN).cls
	time -o $(WORKDIR)/accel_ttime.txt svm_accelerate $(SVMMODEL) $(TRAIN) $(ACCELERATOR)

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

$(MULTICONTROL):
	print_control -Q 6 $(NCLS) > $(MULTICONTROL)

$(TRAIN).vec $(TEST).cls $(TEST).vec: $(TRAIN).cls

#file conversion:
#steps:
#1. copy to working directory
#2. split into test and training, if applicable
#3. normalize and make class labels go from 0 to Nc-1
$(WORKDIR)/heart.trn.cls $(WORKDIR)/segment.trn.cls: $(WORKDIR)/%.trn.cls: $(DATADIR)/%.dat.txt $(WORKDIR)
	lvq2agf -H $< $(WORKDIR)/$(VER))
	agf_preprocess -f $(FRAC) $(WORKDIR)/$(VER) $(TRAIN) $(TEST)
	agf_preprocess -Un $(TRAIN) $(TRAIN)
	agf_preprocess -Un $(TEST) $(TEST)

$(WORKDIR)/sat.trn.cls $(WORKDIR)/shuttle.trn.cls: $(WORKDIR)/%.cls: $(DATADIR)/%.txt $(WORKDIR)
	cp $(DATADIR)/$(VER).trn.txt $(WORKDIR)/$(VER).trn.txt
	cp $(DATADIR)/$(VER).tst.txt $(WORKDIR)/$(VER)).tst.txt
	lvq2agf -H $(WORKDIR)/$(VER).trn.txt $(TRAIN)
	lvq2agf -H $(WORKDIR)/$(VER).tst.txt $(TEST)
	agf_preprocess -Un $(TRAIN) $(TRAIN)
	agf_preprocess -Un $(TEST) $(TEST)

VEHICLE_FILES=xaa.dat xab.dat xac.dat xad.dat xae.dat xaf.dat xag.dat xah.dat xai.dat
VEHICLE_TEST=xaa.dat
VEHICLE_TRAIN=$(VEHICLE_FILES, $(VEHICLE_TEST)=)

$(WORKDIR)/vehicle.trn.cls: $(addprefix, $(DATADIR)/, $(VEHICLE_FILES)) $(WORKDIR)
	cat $(addprefix, $(DATADIR)/, $(VEHICLE_TRAIN)) > $(WORKDIR)/vehicle/vehicle.trn.txt
	cp $(addprefix, $(DATADIR)/, VEHICLE_TEST) $(WORKDIR)/vehicle/vehicle.tst.txt
	sed 's/opel/0; s/saab/1; s/bus/2; s/van/3' $(WORKDIR)/$(VER).trn.txt
	sed 's/opel/0; s/saab/1; s/bus/2; s/van/3' $(WORKDIR)//$(VER).tst.txt
	lvq2agf -H $< $(WORKDIR)/$(VER).trn.txt $(TRAIN)
	lvq2agf -H $< $(WORKDIR)/$(VER).tst.txt $(TEST)
	agf_preprocess -Un $(TRAIN) $(TRAIN)
	agf_preprocess -Un $(TEST) $(TEST)

$(WORKDIR)/covtype.trn.cls: $(WORKDIR)/%.trn.cls: $(DATADIR)/% $(WORKDIR)
	svm2agf -E 0 $(DATADIR)/$(VER) $(WORKDIR)/$(VER)
	agf_preprocess -Unf $(FRAC) $(WORKDIR)/$(VER) $(TRAIN) $(TEST)

$(WORKDIR): $(WORKDIR)
	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
