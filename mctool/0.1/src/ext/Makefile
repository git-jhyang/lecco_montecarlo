
.PHONY: help clean veryclean

include ./Makefile.inc

help :
	#
	# !!! check fortran compiler in 'Makefile' 
	# 
	# help      : print this help
	# lfort     : make library file '$(LIBNAME).a' for fortran
	# linc      : check '$(LIBNAME).a'
	# clean     : remove created files

clean :
	if [ -e $(LIBNAME).a ]; then rm -f $(LIBNAME).a; fi
	rm -f *.o *.mod 

#=================================================================
# Fortran Compiler

###  gFortran

FC       = gfortran
FCFLAGS  = -fno-range-check 



###  intel Fortran

#FC       = ifort
