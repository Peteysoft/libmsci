BASE=letter

SVMMOD=$(BASE).svmmod
AGFMOD=$(BASE).agf
ACCMOD=$(BASE).accel

TEST=$(BASE).tst
TRAIN=$(BASE).trn

#multi-class model:
MODEL=1

#number of classes:
NCLS=26

#fraction to use for testing:
FRAC=0.25

#pre-processing:
PRE=-n

#options for AGF:
AGFOPT=-W 2 -k 200 -t 1e-5 -s 2000 -r-0.9

#options for SVM:
SVMOPT=-h 0 -c 50

#options for SVM-borders:
ACCOPT=-t 1e-5 -s 5000

all: $(TEST).agf.cls $(TEST).svm.cls $(TEST).acc.cls
	@echo "AGF-borders:"
	@echo
	@cls_comp_stats $(TEST).cls $(TEST).agf
	@echo
	@echo "SVM:"
	@echo
	@cls_comp_stats $(TEST).cls $(TEST).svm
	@echo
	@echo "SVM-borders:"
	@echo
	@cls_comp_stats $(TEST).cls $(TEST).acc

clean:
	rm -f $(SVMMOD) $(AGFMOD) $(ACCMOD)
	rm -f $(TEST).*
	rm -f $(TRAIN).*
	#overlaps with first line:
	rm -f $(BASE).svm*
	rm -f $(BASE).agf*
	rm -f $(BASE).acc*

$(TEST).svm.cls $(TEST).svm.con: $(SVMMOD) $(TEST).vec
	classify_m -Z $(SVMMOD) $(TEST).vec $(TEST).svm > $(TEST).svm.txt

$(TEST).agf.cls $(TEST).agf.con: $(AGFMOD) $(TEST).vec
	classify_m $(AGFMOD) $(TEST).vec $(TEST).agf > $(TEST).agf.txt

$(TEST).acc.cls $(TEST).acc.con: $(ACCMOD) $(TEST).vec
	classify_m $(ACCMOD) $(TEST).vec $(TEST).acc > $(TEST).acc.txt

$(SVMMOD): $(TRAIN).svm
	print_control -Q $(MODEL) $(NCLS) | multi_borders -M -- "svm-train -b 1" -+ "$(SVMOPT)" $(TRAIN).svm $(BASE).svm $(SVMMOD)

$(AGFMOD): $(TRAIN).vec $(TRAIN).cls
	print_control -Q $(MODEL) $(NCLS) | multi_borders $(AGFOPT) $(TRAIN) $(BASE).agf $(AGFMOD)

$(ACCMOD): $(TRAIN).vec $(TRAIN).cls $(SVMMOD)
	multi_borders -Z $(ACCOPT) $(SVMMOD) $(TRAIN) $(BASE).acc $(ACCMOD)

$(TRAIN).svm: $(TRAIN).vec $(TRAIN).cls
	agf2ascii -M $(TRAIN) > $(TRAIN).svm

$(TRAIN).vec $(TRAIN).cls: $(BASE).vec $(BASE).cls
	agf_preprocess -zf $(FRAC) $(PRE) letter $(TRAIN) $(TEST)

$(BASE).vec $(BASE).cls: $(BASE).txt
	lvq2agf -H $(BASE).txt $(BASE)

letter.txt: ../statlog/data/letter/letter-recognition.data.txt
	sed 's/A/0/;s/B/1/;s/C/2/;s/D/3/;s/E/4/;s/F/5/;s/G/6/;s/H/7/;s/I/8/;s/J/9/;s/K/10/;s/L/11/;s/M/12/;s/N/13/;s/O/14/;s/P/15/;s/Q/16/;s/Q/16/;s/R/17/;s/S/18/;s/T/19/;s/U/20/;s/V/21/;s/W/22/;s/X/23/;s/Y/24/;s/Z/25/;' < $< | awk -F ',' '{print $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9, $$10, $$11, $$12, $$13, $$14, $$15, $$16, $$1}' > letter.txt


