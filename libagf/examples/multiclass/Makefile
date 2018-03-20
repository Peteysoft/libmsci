SHELL=/bin/bash
TIME=/bin/time

DATAPATH=../datasets

VER=sat

#fraction of data to use for training:
FRAC=0.3
#number of trials:
NTRIAL=10

#number of classes in each dataset:
NCLS_shuttle=7
NCLS_covtype=7
NCLS_covtype_s=7
NCLS_segment=7
NCLS_sat=6
NCLS_poker=10
NCLS_mnist_s=10
NCLS_letter=26
NCLS_pendigits=10
NCLS_usps=10
#this one actually sets the number of classes:
NCLS_humidity=8
NCLS_vehicle=4

NCLS=$(NCLS_$(VER))

#because the data flow is different depending on these parameters, we stick
#all the file naming conventions close together

#type of control file:
#TYPE=1v1		#1 vs. 1
#TYPE=hier		#hierarchical
#TYPE=emp		#empirical hierarchical

#solution method:
#SOL=inv		#matrix inverse
#SOL=iter		#iterative
#SOL=rec		#recursive

WORKDIR=$(VER)
BASE=$(WORKDIR)/$(VER)

KEY = $(TYPE)_$(SOL)
RESULT =$(BASE)_$(KEY)

TRAIN=$(BASE).trn
TEST=$(BASE).tst

#control files:
#training:
TRAINCONTROL = $(BASE)_$(TYPE).mbt

#base name for model files:
MODEL=$(BASE)_$(TYPE)

#classification control file:
#for 1 vs. 1, the solution method is just an option at the classification stage;
#everything else is the same:
CLASSCONTROL_1v1_inv = $(BASE)_1v1
CLASSCONTROL_1v1_iter = $(BASE)_1v1
#for hierarchical classification there are different control files depending
#on the solution method
# there is no need for retraining but the classifiation control file is
# different
# ** marked and un-marked versions... **
CLASSCONTROL_hier_iter = $(BASE)_hier_iter
CLASSCONTROL_hier_rec = $(BASE)_hier
CLASSCONTROL_emp_iter = $(BASE)_emp_iter
CLASSCONTROL_emp_rec = $(BASE)_emp

CLASSCONTROL=$(CLASSCONTROL_$(KEY))

#argument for initial control file:
CONTROL_ARG_1v1 = $(NCLS)
CONTROL_ARG_hier = $(NCLS)
CONTROL_ARG_emp = $(TRAIN)

CONTROL_ARG=$(CONTROL_ARG_$(TYPE))

#command for revised control file:
REVISE_CONTROL_hier_iter = print_control -Q 9 $(BASE)_hier.mbc > $(CLASSCONTROL).mbc
REVISE_CONTROL_emp_iter = print_control -Q 9 $(BASE)_emp.mbc > $(CLASSCONTROL).mbc

REVISE_CONTROL = $(REVISE_CONTROL_$(KEY))

#options for different datasets:
SVMOPT_sat=-c 50 -g 0.1
SVMOPT_segment=-c 100 -g 0.1 
SVMOPT_pendigits=-c 50 -g 0.01
SVMOPT_letter=-c 50 -h 0
SVMOPT_humidity=-h 0 -c 50

OPTIONS=$(SVMOPT_$(VER))

#file extensions for results:
EXT_agf=.vec
EXT_svm=.svmout
EXT=$(EXT_$(METHOD))

#control file option:
TRAINCODE_1v1=5
TRAINCODE_hier=0
TRAINCODE_emp=6

TRAINCODE=$(TRAINCODE_$(TYPE))

#solution option:
TESTCODEinv=5
TESTCODEiter=6
TESTCODErec=0

TESTCODE=$(TESTCODE$(SOL))

#preprocessing options:
PRE=-n

#need this to normalize test data same as training data, if applicable:
ifneq ($(strip $(PRE)),)
  NORMFILE=$(TRAIN).std
  NORMCLAUSE=-a $(NORMFILE)
endif

#subsampling option:
CFLAG_shuttle=correction
#keep class sizes constant:
SAMPLING_METHOD_= -C
#correct for smaller classes:
SAMPLING_METHOD_correction=
#there must be a way of doing it without using three levels of indirection
#(and without resorting to 'if'):
SAMPLING_METHOD=$(SAMPLING_METHOD_$(CFLAG_$(VER)))

#command for subsampling:
#takes into account relative class sizes:
SUBSAMPLE=subsample_special
#doesn't:
#SUBSAMPLE=agf_preprocess

#FCONV_svm=svmout2agf $(RESULT).svmout $(RESULT) > $(RESULT).txt
FCONV=$(FCONV_$(METHOD))

#collect statistics:
STATFILE = $(VER).txt


all:
	make $(BASE)_1v1_inv.txt TYPE=1v1 SOL=inv
	make $(BASE)_1v1_iter.txt TYPE=1v1 SOL=iter
	make $(BASE)_hier_rec.txt TYPE=hier SOL=rec
	make $(BASE)_hier_iter.txt TYPE=hier SOL=iter
	make $(BASE)_emp_rec.txt TYPE=emp SOL=rec
	make $(BASE)_emp_iter.txt TYPE=emp SOL=iter
	cls_comp_stats -Hb $(TEST).cls $(BASE)_1v1_inv | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	cls_comp_stats -Hb $(TEST).cls $(BASE)_1v1_iter | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	cls_comp_stats -Hb $(TEST).cls $(BASE)_hier_rec | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	cls_comp_stats -Hb $(TEST).cls $(BASE)_hier_iter | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	cls_comp_stats -Hb $(TEST).cls $(BASE)_emp_rec | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	cls_comp_stats -Hb $(TEST).cls $(BASE)_emp_iter | tr -d '\n' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	validate_probabilities -H $(TEST).cls $(BASE)_1v1_inv.txt | grep -oE "[+-]?[0-9]+\.[0-9]+(e[+-]?[0-9]+)?" | tr '\n' ' ' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	validate_probabilities -H $(TEST).cls $(BASE)_1v1_iter.txt | grep -oE "[+-]?[0-9]+\.[0-9]+(e[+-]?[0-9]+)?" | tr '\n' ' ' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	validate_probabilities -H $(TEST).cls $(BASE)_hier_iter.txt | grep -oE "[+-]?[0-9]+\.[0-9]+(e[+-]?[0-9]+)?" | tr '\n' ' ' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	validate_probabilities -H $(TEST).cls $(BASE)_emp_iter.txt | grep -oE "[+-]?[0-9]+\.[0-9]+(e[+-]?[0-9]+)?" | tr '\n' ' ' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	validate_probabilities -b $(TEST).cls $(BASE)_hier_rec | grep -oE "[+-]?[0-9]+\.[0-9]+(e[+-]?[0-9]+)?" | tr '\n' ' ' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	validate_probabilities -b $(TEST).cls $(BASE)_hier_iter | grep -oE "[+-]?[0-9]+\.[0-9]+(e[+-]?[0-9]+)?" | tr '\n' ' ' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	validate_probabilities -b $(TEST).cls $(BASE)_emp_rec | grep -oE "[+-]?[0-9]+\.[0-9]+(e[+-]?[0-9]+)?" | tr '\n' ' ' >> $(STATFILE)
	echo -n " " >> $(STATFILE)
	validate_probabilities -b $(TEST).cls $(BASE)_emp_iter | grep -oE "[+-]?[0-9]+\.[0-9]+(e[+-]?[0-9]+)?" | tr '\n' ' ' >> $(STATFILE)
	echo >> $(STATFILE)

#how many lines does results file contain so far?
NLINE!=cat $(STATFILE) | wc -l
stats:
	for ((I=$(NLINE); I<$(NTRIAL); I++)); do \
		make clean; \
		make; \
	done

clean:
	rm -f $(WORKDIR)/*

#performing classifications:
$(RESULT).txt: $(TEST).svm $(CLASSCONTROL).mbaio
	classify_m -Z -Q $(TESTCODE) $(CLASSCONTROL).mbaio $(TEST).vec $(RESULT) > $(RESULT).txt
	$(FCONV)

#compiling the "all-in-one" model file:
$(CLASSCONTROL).mbaio: $(MODEL).mbc
	$(REVISE_CONTROL)
	mbh2mbm -Z $(CLASSCONTROL).mbc $(CLASSCONTROL).mbaio


#training the model:
$(MODEL).mbc: $(TRAIN).svm $(TRAINCONTROL)
	multi_borders -M -- svm-train -+ "-b 1 $(OPTIONS)" \
			$(TRAINCONTROL) $(TRAIN).svm $(MODEL) $(MODEL).mbc

#create control files:
$(TRAINCONTROL): $(TRAIN).svm
	print_control -Q $(TRAINCODE) $(CONTROL_ARG) > $(TRAINCONTROL)

#normalize:
$(TRAIN).svm: $(TRAIN).0.vec
	agf_preprocess $(PRE) $(TRAIN).0 $(TRAIN)
	agf2ascii -M $(TRAIN) > $(TRAIN).svm

$(TEST).svm: $(TEST).0.vec $(TRAIN).svm
	agf_preprocess $(NORMCLAUSE) $(TEST).0 $(TEST)
	agf2ascii -M $(TEST) > $(TEST).svm

#divide into test and training:
$(TRAIN).0.vec $(TEST).0.vec: $(DATAPATH)/$(VER).vec $(WORKDIR)/timestamp
	$(SUBSAMPLE) $(SAMPLING_METHOD) -zf $(FRAC) $(DATAPATH)/$(VER) $(TRAIN).0 $(TEST).0

#make the working directory:
$(WORKDIR)/timestamp:
	if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	date > $(WORKDIR)/timestamp

$(DATAPATH)/$(VER).vec:
	make -C $(DATAPATH)

#$(DATAPATH)/mnist_s.vec:
#	make -C $(DATAPATH) VER=mnist
#	$(SUBSAMPLE) -zf 0.1 $(DATAPATH)/mnist dum $(DATAPATH)/mnist_s

