SBASE=blinking_vortex
RESULT=out.txt
INITCOND=ico.txt
NT=1000
NGRID=100
TSTEP=0.1
SYMBOL=_Z8bv_fielddPfS_Pv

CFLAGS=-shared -fPIC

CC=g++

bv_traj.ps: $(RESULT)
	echo 'set terminal postscript color; set output "bv_traj.ps"; plot "out.txt" using 2:3 with lines' | gnuplot 

$(RESULT): $(INITCOND) $(SBASE).so
	ctraj_model -V 2 -N $(NT) -h $(TSTEP) $(SBASE).so $(SYMBOL) $(INITCOND) $(RESULT) > model_output.txt

tracer.vec tracer.txt: mat.dat
	zonally_symmetric_tracer -V 2 -r 3 | lla2aeb.1 -V 2 -r 3 t0.vec
	echo "tracer.vec=sparse_array(mat.dat)#full(t0.vec)[0]" | sparse_calc

#	echo "sparse_array mat.dat" > sc_script.txt
#	echo "n=size(mat.dat)[1]" >> sc_script.txt
#	echo "t0.vec=empty(n)">>sc_script.txt
#	echo "t0.vec[0:(n/2-1)]=empty(n/2)+1" >> sc_script.txt
#	echo "tracer.vec=mat.dat#t0.vec">>sc_script.txt
#	echo 'print "tracer.txt" tracer.vec'>>sc_script.txt
#	sparse_calc < sc_script.txt

mat.dat: $(SBASE).so
	ctraj_tracer -n ${NGRID} -V 2 -N $(NT) -h $(TSTEP) -r 3 $(SBASE).so $(SYMBOL) mat.dat

$(SBASE).so: $(SBASE).cc
	$(CC) $(CFLAGS) $(SBASE).cc -o $(SBASE).so

