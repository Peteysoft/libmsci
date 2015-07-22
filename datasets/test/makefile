VPATH = ../

CC = g++
CFLAGS = -g -Wno-deprecated

test_interpol: test_interpol.cc libdataset-g.a
	$(CC) $(CFLAGS) -o test_interpol test_interpol.cc ../libdataset-g.a \
		-I ../../include

test_interpol_old: test_interpol_old.cc libdataset-g.a ../old/interpol_old.cc
	$(CC) $(CFLAGS) -o test_interpol_old ../old/interpol_old.cc \
		test_interpol_old.cc ../libdataset-g.a \
		-I ../../include

dep.exe: dependent_interface.cpp dataset_lib.a dependent_temp.h
	$(CC) $(CFLAGS) -o dep.exe dependent_interface.cpp ../dataset_lib.a

simple.exe: simple_interface.cpp simple_temp.h dataset_lib.a
	$(CC) $(CFLAGS) -o simple.exe simple_interface.cpp ../dataset_lib.a

dataset_lib.a:
	make -C ../ dataset_lib.a

