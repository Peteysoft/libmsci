#EXE_EXT = .exe

#DATATYPE=ncep
DATATYPE=badc
#DATATYPE=ecmwf

#command-line calculator:
CCALC = bc

#DATAPATH = /cygdrive/e/data/ncep/sixhour
#DATAPATH = /media/Lexar/ncep_data
DATAPATH = /home/peter/data/$(DATATYPE)_data

BINPATH = /home/peter/bin

# so we can have the files in a different working directory:
WORKPATH = .
#VFIELD_PATH = /cygdrive/e/tmpdata/vfield/
VFIELD_PATH = .
#VFIELD_PATH = /tmp/vfield/
#WORKPATH = /work/patmills/ctraj_work
YEAR = 1999
NEXT_YEAR=$(shell echo $(YEAR)+1 | $(CCALC))
LAST_YEAR=$(shell echo $(YEAR)-1 | $(CCALC))
LEV = 530

#date range of NCEP velocity data to import
T1= $(YEAR)/1/1
T2= $(YEAR)/2/31

T0 = $(T1)

COPTS = -z 100 -I 0 -F 1

BASE = $(YEAR)_$(LEV)K
TFLAG = -T

#resolution of wind fields
VGRID = 200

#for tracer simulation:
NT = 20		#number of time steps
NT2 = $(NT)
NGRID=100	#resolution of tracer field
TBASE = $(WORKPATH)/t$(BASE)

#for contour advection simulation:
ARC = 1		#maximum degrees of arc
MINS = 10	#minimum spacing between pairs of points
MAXS = 100	#maximum spacing between pairs of points
NC = 20		#number of time steps
C0LAT = 45
CBASE = $(WORKPATH)/c$(BASE)_$(C0LAT)deg

VFIELD_N = $(VFIELD_PATH)/vfield_$(DATATYPE)$(BASE)_N.ds
VFIELD_S = $(VFIELD_PATH)/vfield_$(DATATYPE)$(BASE)_S.ds


VPATH = $(BINPATH) $(DATAPATH) $(WORKPATH) ../plotting

DATA_FILES = $(DATAPATH)/uwnd.$(YEAR).nc $(DATAPATH)/vwnd.$(YEAR).nc $(DATAPATH)/air.$(YEAR).nc \
	$(DATAPATH)/uwnd.$(NEXT_YEAR).nc $(DATAPATH)/vwnd.$(NEXT_YEAR).nc $(DATAPATH)/air.$(NEXT_YEAR).nc \
	$(DATAPATH)/uwnd.$(LAST_YEAR).nc $(DATAPATH)/vwnd.$(LAST_YEAR).nc $(DATAPATH)/air.$(LAST_YEAR).nc

TINIT = $(TBASE).init
MATFILE = $(WORKPATH)/mat$(BASE).dat
TFILE = $(TBASE).vec

CINIT = $(CBASE).init
C0FILE = $(WORKPATH)/c$(C0LAT)deg.bev
CFILE = $(CBASE).bev

TANIM = $(TBASE).gif
CANIM = $(CBASE).gif

#commands:
TRACER = ctraj_tracer
CONTOUR = ctraj_contour

both: $(CFILE) $(TFILE) animate_combined.sh gen_zgrid
	animate_combined.sh $(COPTS) $(TFILE) $(CFILE) $(WORKPATH)/b$(BASE).gif

contour: $(CANIM)

tracer: $(TANIM)

vfield: $(VFIELD_N) $(VFIELD_S)

clean:
	rm -f $(CANIM)
	rm -f $(TANIM)
	rm -f $(TINIT)
	rm -f $(MATFILE)
	rm -f $(TFILE)
	rm -f $(CINIT)
	rm -f $(C0FILE)
	rm -f $(CFILE)
	rm -f zsym.dat

data:
	touch $(DATA_FILES)
	make $(VFIELD_N) $(VFIELD_S)

clean_all:
	make clean
	rm -f $(VIELD_N)
	rm -f $(VIELD_S)

##############################################################
# this section for building the contour advection simulation
# ############################################################

$(CANIM): $(CFILE) animate_ctraj.sh
	animate_ctraj.sh $(CFILE) $(CANIM)

$(CFILE): $(CONTOUR)$(EXE_EXT) $(CINIT) $(VFIELD_N) $(VFIELD_S) $(C0FILE)
	#contour2 $(CINIT)
	#cp $(C0FILE) $(CFILE)
	$(CONTOUR) -c $(ARC) -m $(MINS) -s $(MAXS) -i $(T0) -N $(NC) -h 1 -k 4 \
		$(VFIELD_N) $(VFIELD_S) $(C0FILE) $(CFILE)

$(CINIT): Makefile
	echo $(VFIELD_N) > $@
	echo $(VFIELD_S) >> $@
	echo $(C0FILE) >> $@
	echo $(CFILE) >>$@
	echo $(T0) >> $@
	echo 1 >> $@
	echo $(NC) >> $@
	echo 4 >> $@
	echo $(ARC) >> $@
	echo $(MINS) >> $@
	echo $(MAXS) >> $@

$(C0FILE): c0$(EXE_EXT) xy2bev$(EXE_EXT)
	c0 $(C0LAT) | xy2bev $(C0FILE)

##############################################################
# this section for building the tracer simulation
# ############################################################

# animate the field:
$(TANIM): $(TFILE) animate_ctraj.sh 
	animate_ctraj.sh $(COPTS) -q -N $(NT2) $(TFILE) $(TANIM)

# calculate actual tracer fields:
$(TFILE): $(MATFILE) sparse_vect_prod$(EXE_EXT) zsym.dat
	cp zsym.dat $(TFILE)
	sparse_vect_prod -n $(NT) $(MATFILE) $(TFILE)

# create matrix mapping:
$(MATFILE): $(TINIT) $(TRACER)$(EXE_EXT) $(VFIELD_N) $(VFIELD_S)
	#tracer4_matout $(TINIT)
	$(TRACER) -n $(NGRID) -h 1 -k 4 -i $(T0) -N $(NT) $(VFIELD_N) $(VFIELD_S) $(MATFILE)

# create initialization file:
$(TINIT): Makefile
	echo $(VFIELD_N) > $@
	echo $(VFIELD_S) >> $@
	echo $(MATFILE) >> $@
	echo mapout.txt >>$@
	echo $(T0) >> $@
	echo 1 >> $@
	echo $(NT) >> $@
	echo 4 >> $@
	echo $(NGRID) >> $@

# create zonally-symmetric tracer:
zsym.dat: zonally_symmetric_tracer$(EXE_EXT) lla2aeb$(EXE_EXT) Makefile
	zonally_symmetric_tracer | lla2aeb -n $(NGRID) zsym.dat

################################################
# this section for converting the velocity fields
# ##############################################

$(VFIELD_PATH)/vfield_ncep$(BASE)_N.ds $(VFIELD_PATH)/vfield_ncep$(BASE)_S.ds: $(DATA_FILES) #nc2ds$(EXE_EXT)
	nc2ds -n $(VGRID) $(TFLAG) -p $(DATAPATH) \
		-i $(LAST_YEAR)/12/31 -f $(NEXT_YEAR)/1/1 $(LEV) \
		$(VFIELD_S) $(VFIELD_N)

#$(VFIELD_S): $(DATA_FILES) #nc2ds$(EXE_EXT)
#	nc2ds -n $(VGRID) $(TFLAG) -- -p $(DATAPATH) \
#		-i $(LAST_YEAR)/12/31 -f $(NEXT_YEAR)/1/1 $(LEV) $(VFIELD_S)

$(VFIELD_PATH)/vfield_ecmwf$(BASE)_N.ds $(VFIELD_PATH)/vfield_ecmwf$(BASE)_S.ds: $(DATA_FILES) #ecmwf2ds.sh
	ecmwf2ds -n $(VGRID) $(TFLAG) -p $(DATAPATH) -O \
		-i $(LAST_YEAR)/12/31 -f $(NEXT_YEAR)/1/1 $(LEV) \
		$(VFIELD_S) $(VFIELD_N)

$(VFIELD_PATH)/vfield_badc$(BASE)_N.ds $(VFIELD_PATH)/vfield_badc$(BASE)_S.ds:
	import_vfield.sh -V 0 -n $(VGRID) $(TFLAG) \
		-i $(LAST_YEAR)/12/31-12 -f $(NEXT_YEAR)/1/1-12 $(LEV) \
		$(VFIELD_S) $(VFIELD_N)

