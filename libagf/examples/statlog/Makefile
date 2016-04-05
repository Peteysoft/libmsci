V=2

VER1=australian
VER2=german

NCLS1=2
NCLS2=2
NCLS=$(NCLS$(V))

AGFOPT1=
AGFOPT2=
AGFOPT=$(AGFOPT$(V))

SCRIPTDIR=../../scripts

DATADIR=statlog_data/$(VER$(V))

WORKDIR=work/$VER$(V)

TRAINFILE0=$(DATADIR)/$(VER$(V))

BASE=$(WORKDIR)/$(VER$(V))
TRAINFILE=$(BASE)

MULTICONTROL=onevsone.mbc

SVMMODEL=$(BASE).svmmod

$(WORKDIR)/nbv.txt: $(MULTICONTROL) $(TRAINFILE)
	$(SCRIPTDIR)/mb_nv.sh $(AGFOPT) $(MULTICONTROL) $(TRAINFILE) $(WORKDIR)/nbv.txt

$(SVMMODEL): $(MULTICONTROL) $(TRAINFILE)
	agf2ascii -M $(TRAINFILE) $(TRAINFILE).svm
	multi_borders -- "svm-multi $(SVMOPT)" $(MULTICONTROL) $(TRAINFILE) $(SVMMODEL) &> svmmulti.log

$(TRAINFILE): $(TRAINFILE0)
	lvq2agf $(TRAINFILE0) $(TRAINFILE)

$(MULTICONTROL):
	print_control -Q 6 $(NCLS) > $(MULTICONTROL)

all:
	echo $(DIR)
