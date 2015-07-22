#EXE_EXT = .exe

BIN_PATH=/usr/local/bin

all: gen_zgrid$(EXE_EXT)
	install gen_zgrid plot_frame.sh animate_ctraj.sh animate_combined.sh \
			quick_plot.sh $(BIN_PATH)

gen_zgrid$(EXE_EXT): gen_zgrid.cc
	g++ gen_zgrid.cc -o gen_zgrid
