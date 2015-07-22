T0=2000/1/1
TF=2000/4/1

NGRID=50
LEV=500
YEAR=2000

VFIELD_PATH = /home/petey/data/vfield
CTRAJ_TEST_PATH=/home/petey/peteysoft/ctraj/scripts/test

BASE=$(YEAR)_$(LEV)K
VFIELD_N = $(VFIELD_PATH)/vfield_ncep$(BASE)_N.ds
VFIELD_S = $(VFIELD_PATH)/vfield_ncep$(BASE)_S.ds

tforward.dat: $(VIELD_N) $(VFIELD_S)
	ctraj_tracer -n $(NGRID) -i $(T0) -f $(TF) $(VFIELD_S) $(VFIELD_N) \
		tforward.dat

tbackward.dat: $(VIELD_N) $(VFIELD_S)
	ctraj_tracer -n $(NGRID) -i $(TF) -f $(T0) $(VFIELD_S) $(VFIELD_N) \
		tbackward.dat

$(VFIELD_N) $(VFIELD_S):
		make -C $(CTRAJ_TEST_PATH) vfield

