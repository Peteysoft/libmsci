#!/bin/bash
cd /home/home/pmills/my_software/ctraj/scripts
make -f Makefile.proxy clean
(time make -f Makefile.proxy) &> proxymake_err.txt

