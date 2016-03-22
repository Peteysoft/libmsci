#drive letter containing Landsat data
DRIVELETTER=e
NVEC=63			#number of features in (raw) training data

NSV=10			#number of singular vectors to keep

#parameters for classication scheme
#CLASS_OPTS=-c 1 -W 3		#direct classification

CLASS_OPTS=-c 5 -W 3 control2.mbt	#multi-borders classification

#path for python scripts:
#BINPATH=/home/Petey/my_programs/libagf/applications/Landsat
BINPATH=.

#path for l1 Landsat raster (tiff) data
DATAPATH=/cygdrive/$(DRIVELETTER)/data/Landsat
#DATAPATH=/media/petey/1A8F-0227/Landsat

TRAINING0=clearcut_class_data
TRAINING1=cccd1

#nfold cross-validation for parameter optimization:
test: $(TRAINING1).vec $(TRAINING1).cls
	nfold $(CLASS_OPTS) $(TRAINING1)

clean:
	rm -f $(TRAINING1).vec $(TRAINING1).cls $(TRAINING1).std
	rm -f $(TRAINING0).vec $(TRAINING0).cls $(TRAINING0).lvq
	rm -f *_clearcut_t*.lvq

#pre-processing: normalization, SVD, features pruning
$(TRAINING1).vec $(TRAINING1).cls: $(TRAINING0).vec $(TRAINING0).cls Makefile
	agf_preprocess -S $(NSV) $(TRAINING0) $(TRAINING1)

#convert to binary file format:
$(TRAINING0).vec $(TRAINING0).cls: $(TRAINING0).lvq
	lvq2agf $(TRAINING0).lvq $(TRAINING0)

#extract features data from Landsat scenes and compile into a single
#ASCII file:
$(TRAINING0).lvq:
	echo $(NVEC) > $(TRAINING0).lvq
	make LVQFile LSBASE=LT50190272011278EDC00 TRIAL=1
	make LVQFile LSBASE=LT50190272011278EDC00 TRIAL=2
	make LVQFile LSBASE=LT50190272011278EDC00 TRIAL=3
	cat *_clearcut_t*.lvq >> $(TRAINING0).lvq

LVQFile: $(LSBASE)_clearcut_t$(TRIAL).lvq

#extract features data from image coordinates:
$(LSBASE)_clearcut_t$(TRIAL).lvq:  $(LSBASE)_clearcut_t$(TRIAL).txt
	$(BINPATH)/extract_landsat_features.py --texture $(DATAPATH)/$(LSBASE)/$(LSBASE) $< > $@ 

