
AR = ar

# Intel compiler needed the following definition for the archiver:
# AR = xiar

ARFLAGS = r

RANLIB = ranlib

# include norm2.o below, if your compiler is missing that intrinsic
OBJECTS = ReadHeader.o beampatternMod.o MathConstantsMod.o RefCoMod.o SdRdRMod.o sspMod.o WriteHeader.o errout.o interpMod.o MergeVectors.o munk.o ReadEnvironment.o SortMod.o splinec.o subtabulate.o twersky.o weight.o xerror.o norms.o cross_products.o monotonicMod.o PolyMod.o crci.o PekerisRoot.o

all: libmisc.a

install: all

clean:
	-rm -f *.a *.mod *.o *_genmod.f90

# Declare implicit rules (for make utilities with missing or broken build rules)

.SUFFIXES: .a .o .mod .f .f90

.f.o:
	$(FC) -c $(FFLAGS) $< 

.f90.o:
	$(FC) -c $(FFLAGS) $< 

.f90.mod:
	$(FC) -c $(FFLAGS) $< 


# Declare explicit dependencies (and explicit build rules as needed)

libmisc.a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)
	$(RANLIB) $@
	@echo "Routines needed by SCOOTER\SPARC\BELLHOP\BOUNCE built"

subtabultate.o subtabulate.mod: subtabulate.f90
	$(FC) -c $(FFLAGS) $<

SdRdRMod.o sdrdrmod.mod: SdRdRMod.f90 monotonicmod.mod sortmod.mod subtabulate.mod
	$(FC) -c $(FFLAGS) $< 

ReadHeader.o: sdrdrmod.mod ReadHeader.f90

PolyMod.o polymod.mod: PolyMod.f90
	$(FC) -c $(FFLAGS) $< 

crci.o: crci.f90 MathConstantsMod.mod
	$(FC) -c $(FFLAGS) $<

beampatternMod.o beampatternmod.mod: beampatternMod.f90
	$(FC) -c $(FFLAGS) $<

interpMod.o interpmod.mod: interpMod.f90
	$(FC) -c $(FFLAGS) $<

MathConstantsMod.o mathconstantsmod.mod: MathConstantsMod.f90
	$(FC) -c $(FFLAGS) $< 

RefCoMod.o refcomod.mod: RefCoMod.f90 polymod.mod
	$(FC) -c $(FFLAGS) $< 

sspMod.o sspmod.mod: sspMod.f90
	$(FC) -c $(FFLAGS) $< 

monotonicMod.o monotonicmod.mod: monotonicMod.f90
	$(FC) -c $(FFLAGS) $< 

SortMod.o sortmod.mod: SortMod.f90
	$(FC) -c $(FFLAGS) $< 

PekerisRoot.o: PekerisRoot.f90

WriteHeader.o: WriteHeader.f90 SdRdRMod.mod
errout.o: errout.f90
MergeVectors.o: MergeVectors.f90
munk.o: munk.f90
norms.o norms.mod: norms.f90
	$(FC) -c $(FFLAGS) $< 


ReadEnvironment.o: ReadEnvironment.f90 sspMod.mod
splinec.o: splinec.f90
twersky.o: twersky.f90
weight.o: weight.f90
xerror.o: xerror.f90
cross_products.o cross_products.mod: cross_products.f90
	$(FC) -c $(FFLAGS) $< 

