MINSAMPLE=5
MAXSAMPLE=200
NTEST=10
NTRIAL=10

V=1

VER1=australian
VER2=german
VER3=heart

NCLS1=2
NCLS2=2
NCLS3=2

NCLS=$(NCLS$(V))

AGFOPT1=
AGFOPT2="-W 2"
AGFOPT3=

AGFOPT=$(AGFOPT$(V))

SVMOPT1=
SVMOPT2=
SVMOPT3=

SVMOPT=$(SVMOPT$(V))

SCRIPTDIR=../../scripts

DATADIR=statlog_data/$(VER$(V))

WORKDIR=work/VER$(V)

TRAINFILE0=$(DATADIR)/$(VER$(V)).dat.txt

BASE=$(WORKDIR)/$(VER$(V))
TRAIN=$(BASE)

MULTICONTROL=onevsone.mbc

SVMMODEL=$(BASE).svmmod
AGFMODEL=$(BASE).agf
ACCELERATOR=$(BASE).accel

NBORDER_AGF=$(WORKDIR)/nborder_agf.txt
NBORDER_SVM=$(WORKDIR)/nborder_svm.txt
NSUPPORT=$(WORKDIR)/nsupport.txt

all: $(NBORDER_SVM) $(NBORDER_AGF) $(NSUPPORT)

clean:
	rm -f $(BASE).*

#find number of border samples SVM-borders:
$(NBORDER_SVM): $(SVMMODEL).svc $(TRAIN).vec $(TRAIN).cls
	$(SCRIPTDIR)/mb_nb.sh -g -q $(NTEST) -N $(NTRIAL) -M \
		-s $(MINSAMPLE) -S $(MAXSAMPLE) \
		$(SVMMODEL).svmmod $(TRAIN) $(NBORDER_SVM)

#find number of border samples AGF-borders:
$(NBORDER_AGF): $(MULTICONTROL) $(TRAIN).vec $(TRAIN).cls
	$(SCRIPTDIR)/mb_nb.sh -g -q $(NTEST) -N $(NTRIAL) \
		-s $(MINSAMPLE) -S $(MAXSAMPLE) $(AGFOPT) \
		$(MULTICONTROL) $(TRAIN) $(NBORDER_AGF)

$(NSUPPORT): $(TRAIN).svm
	$(SCRIPTDIR)/svm_nsv.sh -G -q $(NTEST) $(SVMOPT) $(TRAIN) $(NSUPPORT)

$(ACCELERATOR): $(MULTICONTROL) $(SVMMODEL).svc $(TRAIN).vec $(TRAIN).cls
	multi_borders -Z $(SVMMODEL).svc $(TRAIN) $(ACCELERATOR) $(ACCELERATOR).mbc &>> svmmulti.log

$(SVMMODEL).svmmod: $(TRAIN).svm
	svm-predict -b 1 $(SVMOPT) $(TRAIN).svm $(SVMMODEL).svmmod

$(SVMMODEL).svc: $(MULTICONTROL) $(TRAIN).svm
	multi_borders -M -- "svm-train $(SVMOPT)" $(MULTICONTROL) $(TRAIN).svm $(SVMMODEL) $(SVMMODEL).svc &> svmmulti.log

$(TRAIN).svm: $(TRAIN).cls $(TRAIN).vec
	agf2ascii -M $(TRAIN) $(TRAIN).svm

$(TRAIN).cls $(TRAIN).vec: $(TRAINFILE0)
	lvq2agf -H $(TRAINFILE0) $(TRAIN)
	agf_preprocess -Un $(TRAIN) $(TRAIN)

$(MULTICONTROL):
	print_control -Q 6 $(NCLS) > $(MULTICONTROL)

