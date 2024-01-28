
.PHONY: all lfort linc
.SUFFIXES: .f90 .o


LIBNAME  = liblocfort

#=================================================================

SRC      = constants util iotool vasp
SRC_MT   = mt_serial

#=================================================================

OBJS      = $(addsuffix .o, ${SRC})
OBJS_MT   = $(addsuffix .o, ${SRC_MT})

OBJS_ALL  = $(OBJS) $(OBJS_MT)

#=================================================================

constants.o : 
mt_serial.o : 
iotool.o    : constants.o
util.o 	    : constants.o iotool.o
vasp.o      : constants.o iotool.o util.o

#=================================================================

all : lfort

lfort : $(LIBNAME).a 

linc :
	ar t $(LIBNAME).a

#=======================================================================

$(LIBNAME).a : $(OBJS_ALL)
	ar src $@ $(OBJS_ALL)

%.o : %.f90
	$(FC) -c $(FCFLAGS) $< -o $@

