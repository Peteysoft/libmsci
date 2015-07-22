#EXE_EXT = .exe

#DATAPATH = /data/ncep/sixhour
DATAPATH = /home/home/pmills/tracer_correlation
#VFIELD_PATH = /temp/vfield
VFIELD_PATH = .

#BINPATH = /home/petey/bin
BINPATH = /home/home/pmills/bin

# so we can have the files in a different working directory:
WORKPATH = work
PREV = 1998
YEAR = 1999
NEXT = 2000
LEV = 500

#date range of NCEP velocity data to import
T1= $(PREV)/12/31-18
T2= $(NEXT)/1/1

T0 = $(T1)

BASE = $(YEAR)_$(LEV)K
TFLAG = -T

#resolution of wind fields
VGRID = 120

#for contour advection simulation:
DSMIN = 1	#minimum spacing between pairs of points
DSMAX = 3  	#maximum spacing between pairs of points
WRITE_INT=4
NC = 20		#number of write intervals 
I0 = 0
N = 140
MAXNPT = 500000

BS = 16777216

C0LAT = 45
CSTART = $(WORKPATH)/c$(C0LAT)deg.bev
CWORK = $(WORKPATH)/cswap.bev2
CBASE = $(WORKPATH)/cspace

VFIELD_N = $(VFIELD_PATH)/vfield_ncep$(BASE)_N.ds
VFIELD_S = $(VFIELD_PATH)/vfield_ncep$(BASE)_S.ds

VPATH = $(BINPATH) $(DATAPATH) $(VFIELD_PATH)

DATA_FILES = $(DATAPATH)/uwnd.$(YEAR).nc $(DATAPATH)/vwnd.$(YEAR).nc $(DATAPATH)/air.$(YEAR).nc \
	$(DATAPATH)/uwnd.$(NEXT).nc $(DATAPATH)/vwnd.$(NEXT).nc $(DATAPATH)/air.$(NEXT).nc \
	$(DATAPATH)/uwnd.$(PREV).nc $(DATAPATH)/vwnd.$(PREV).nc $(DATAPATH)/air.$(PREV).nc

unc_exp_: 
	./run_unc_exp2.sh -+ -u 1000 -U 100000 -I 0.001 -e 40 -g $(CBASE) $(I0) $(N) $(WRITE_INT)
#	./run_bev2vec.sh $(CBASE)

advect: $(CBASE).bev $(VFIELD_N) $(VFIELD_S)
	./run_contour2.sh -O $(WRITE_INT) -M $(MAXNPT) $(VFIELD_N) $(VFIELD_S) $(CBASE) $(I0) $(N)

advect_old: $(CWORK) $(VFIELD_N) $(VFIELD_S)
	./run_contour3.sh -N $(WRITE_INT) -m $(DSMIN) -s $(DSMAX) $(VFIELD_N) $(VFIELD_S) $(CWORK) $(I0) $(NC) $(CBASE)

space_filling.txt:
	./run_unc_exp.sh -+ -N $(WRITE_INT) $(CBASE) $(I0) $(NC) > space_filling.txt

$(CWORK): $(CSTART)
	bev2bev2 -B $(BS) $(CSTART) $(CWORK)
	rm $(CSTART)

$(CBASE).bev: 
	c0 $(C0LAT) | xy2bev $(CBASE).bev

vfield: $(VFIELD_S) $(VFIELD_N)

################################################
# this section for converting the velocity fields
# ##############################################

$(VFIELD_N): $(DATA_FILES) #nc2ds$(EXE_EXT)
	nc2ds -n $(VGRID) $(TFLAG) -p $(DATAPATH) \
		-i $(T1) -f $(T2) $(LEV) $(VFIELD_N)

$(VFIELD_S): $(DATA_FILES) #nc2ds$(EXE_EXT)
	nc2ds -n $(VGRID) $(TFLAG) -- -p $(DATAPATH) \
		-i $(T1) -f $(T2) $(LEV) $(VFIELD_S)

