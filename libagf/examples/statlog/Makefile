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

NCLS=$(NCLS_$(VER))

#special pre-processing:
PRE_heart=-n
PRE_shuttle=-n
PRE_covtype=-n
PRE_segment=-n
PRE_sat=-nS 9
PRE_poker=-n
PRE_vehicle=-n
#PRE_mnist=-nS 20

PRE=$(PRE_$(VER))

#options for AGF-borders:
AGFOPT_heart=-W 5
AGFOPT_shuttle=-W 1.5 -k 200 -s 250
AGFOPT_covtype=
AGFOPT_segment=-W 3 -s 100 -k 100
AGFOPT_sat=-W 5 -s 200 -k 100
AGFOPT_poker=-k 200 -W 2
AGFOPT_mnist=-k 400 -W 25
AGFOPT_vehicle=-W 5 -s 100 -k 100

AGFOPT=$(AGFOPT_$(VER))

#options for svm-train:
SVMOPT_heart=-c 0.5 -g 0.01
SVMOPT_sat=-c 50 -g 0.1
SVMOPT_segment=-c 100 -g 0.1 
SVMOPT_vehicle=-c 100 -g 0.01

SVMOPT=$(SVMOPT_$(VER))

#number of nearest neighbours for knn:
K_heart=-k 11
K_shuttle=-k 1
K_covtype=-k 11
K_segment=-k 3
K_sat=-k 15
K_vehicle=-k 9

K=$(K_$(VER))

#number of border samples (or other options) for "accelerated" SVM model:
NSAMP_segment=-s 50
NSAMP_vehicle=-s 50
NSAMP_sat=-s 150

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

RESULTS=results.txt

MULTICONTROL=$(WORKDIR)/onevsone.mbc

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

#phony targets:
.PHONY: all
.PHONY: clean
.PHONY: nborders

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

stats: calc_stats
	rm -f $(VER).txt
	for I in {1..$(NRUN)}; do \
		make clean; \
		make RESULTS=$(VER).txt; \
	done
	./calc_stats $(VER).txt

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
	@cls_comp_stats -Hb $(TEST).cls $(TEST).knn | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	@cls_comp_stats -Hb $(TEST).cls $(TEST).agf | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	@cls_comp_stats -Hb $(TEST).cls $(TEST).svm | tr -d '\n' >> $(RESULTS)
	echo -n " " >> $(RESULTS)
	@cls_comp_stats -Hb $(TEST).cls $(TEST).acc >> $(RESULTS)

nborders: $(NBORDER_SVM)

clean:
	rm -f $(WORKDIR)/*

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
$(TEST).knn.cls $(TEST).knn.cor: $(TRAIN).vec $(TRAIN).cls $(TEST).vec
	time -o $(WORKDIR)/knn_time.txt knn $K classify $(TRAIN) $(TEST).vec $(TEST).knn

$(TEST).svm.cls $(TEST).svm.cor: $(SVMMODEL) $(TEST).svm
	time -o $(WORKDIR)/svm_time.txt svm-predict -b 1 $(TEST).svm $(SVMMODEL) $(TEST).svm.txt
	svmout2agf $(TEST).svm.txt $(TEST).svm

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
	print_control -Q 6 $(NCLS) | time -o $(WORKDIR)/agf_ttime.txt multi_borders $(AGFOPT) $(TRAIN) $(AGFMODEL) $(AGFMODEL).mbc

#only used for one script:
$(MULTICONTROL):
	print_control -Q 6 $(NCLS) > $(MULTICONTROL)

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
	agf_preprocess $(PRE) -Uzf $(FRAC) $(WORKDIR)/$(VER) $(TRAIN) $(TEST)

#3rd feature must be removed before normalization 
#since it only takes on one value:
$(WORKDIR)/segment.trn.cls: $(DATADIR)/segment.dat.txt
	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	lvq2agf -H $< $(WORKDIR)/$(VER)
	agf_preprocess $(PRE) -FzUf $(FRAC) $(WORKDIR)/$(VER) $(TRAIN) $(TEST) 0 1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18

$(WORKDIR)/sat.trn.cls $(WORKDIR)/shuttle.trn.cls: $(WORKDIR)/%.cls: $(DATADIR)/%.txt
	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	lvq2agf -H $(DATADIR)/$(VER).trn.txt $(TRAIN)
	lvq2agf -H $(DATADIR)/$(VER).tst.txt $(TEST)
	agf_preprocess -U $(PRE) $(TRAIN) $(TRAIN)
	agf_preprocess $(NORMCLAUSE) -U $(TEST) $(TEST)

VEHICLE_SUFFIX=.dat.txt
VEHICLE_FILES=xaa xab xac xad xae xaf xag xah xai
VEHICLE_TEST=xaa
VEHICLE_TRAIN=$(addsuffix $(VEHICLE_SUFFIX), $(filter-out $(VEHICLE_TEST), $(VEHICLE_FILES)))

#$(WORKDIR)/vehicle.$(VEHICLE_TEST).trn.cls: $(WORKDIR)/vehicle.trn.cls \

$(WORKDIR)/vehicle.trn.cls: $(addprefix, $(DATADIR)/, $(VEHICLE_FILES))
	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	cat $(addprefix $(DATADIR)/, $(VEHICLE_TRAIN)) | sed 's/opel/0/;s/saab/1/;s/bus/2/;s/van/3/' > $(WORKDIR)/vehicle.trn.txt
	cat $(DATADIR)/$(VEHICLE_TEST)$(VEHICLE_SUFFIX) | sed 's/opel/0/;s/saab/1/;s/bus/2/;s/van/3/' > $(WORKDIR)/$(VER).tst.txt
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

