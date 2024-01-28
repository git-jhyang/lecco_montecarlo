
.SUFFIXES: .f90 .o .x
.PHONY: config test opt gen

# target ============================================================

LOC_BIN = ../../bin
TG_TEST = test.x
TG_OPT  = opt.x
TG_GEN  = gen.x

# objects ==========================================================

OBJS_OPT =  parallel.o input.o lattice.o coordinate.o optimizer.o energy.o
OBJS_GEN =  parallel.o input.o coordinate.o generator.o
OBJS_ALL =  $(OBJS_OPT) generator.o

parallel.o	:
input.o		: parallel.o
lattice.o	: parallel.o input.o 
coordinate.o	: parallel.o input.o
energy.o	: parallel.o input.o coordinate.o
optimizer.o	: parallel.o input.o lattice.o coordinate.o energy.o
generator.o	: parallel.o input.o coordinate.o 
main_opt.o	: $(OBJS_OPT)
main_gen.o	: $(OBJS_GEN)
main_test.o	: $(OBJS_ALL)

# libraries =========================================================
# basic extension

LOC_EXT = ./ext
#LOC_MTST = $(LOC_EXT)/mt_stream_f90-1.11

LIB_EXT = -L$(LOC_EXT) -llocfort
INC_EXT = -I$(LOC_EXT) #-I$(LOC_MTST)

# aenet extension

LOC_AENET = ./aenet-1.0.0_src

LIB_AENET = -L$(LOC_AENET) -laenet 
INC_AENET = -I$(LOC_AENET)

# gathering 

ALLLIB = $(LIB_EXT) $(LIB_MKL) $(LIB_AENET)
ALLINC = $(INC_EXT) $(INC_MKL) $(INC_AENET)

# rules ============================================================

%.o : %.f90
	$(MPIFC) -c $(MPIFLAGS) $< $(LIB_EXT) $(INC_EXT) -o $@

energy.o : energy.f90
	$(MPIFC) -c $(MPIFLAGS) $< $(LIB_EXT) $(LIB_AENET) $(INC_EXT) $(INC_AENET) -o $@

$(LOC_BIN)/$(TG_OPT) : main_opt.o
	$(MPIFC) $(MPIFLAGS) $< $(OBJS_OPT) $(ALLLIB) $(ALLINC) -o $@

$(LOC_BIN)/$(TG_GEN) : main_gen.o
	$(MPIFC) $(MPIFLAGS) $< $(OBJS_GEN) $(LIB_EXT) $(LIB_MKL) $(INC_EXT) $(INC_MKL) -o $@

$(LOC_BIN)/$(TG_TEST) : main_test.o
	$(MPIFC) $(MPIFLAGS) $< $(OBJS_ALL) $(ALLLIB) $(ALLINC) -o $@

# ==================================================================

test : $(LOC_BIN)/$(TG_TEST)
opt : $(LOC_BIN)/$(TG_OPT)
gen : $(LOC_BIN)/$(TG_GEN)

all : config opt gen

