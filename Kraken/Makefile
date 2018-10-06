
EXECUTABLES = bounce.exe kraken.exe krakenc.exe

BOUNCE_OBJECTS = bounce.o BCImpedancec.o KrakencMod.o kuping.o ../misc/libmisc.a
BOUNCE_LIBS = -L../misc -lmisc

KRAKEN_OBJECTS = BCImpedance.o KrakenMod.o kuping.o InverseIterationMod.o RootFinderSecantMod.o zbrentx.o kraken.o ../misc/libmisc.a
KRAKEN_LIBS = -L../misc -lmisc

KRAKENC_OBJECTS = BCImpedancec.o KrakencMod.o kuping.o InverseIterationMod.o RootFinderSecantMod.o krakenc.o ../misc/libmisc.a
KRAKENC_LIBS = -L../misc -lmisc

all: $(EXECUTABLES)
	@echo "KRAKEN built"
	@echo "************"
	@echo " "

install: all
	for f in $(EXECUTABLES) ; do \
		echo "----- Installing $$f"; cp -p $$f ../bin; \
	done

clean:
	-rm -f *.o *.mod *.exe *_genmod.f90

# Declare implicit rules (for make utilities with missing or broken build rules)

.SUFFIXES: .a .o .mod .f .f90

.f.o:
	$(FC) -c $(FFLAGS) $< 

.f90.o:
	$(FC) -c $(FFLAGS) $< 

.f90.mod:
	$(FC) -c $(FFLAGS) $< 

# Declare explicit dependencies (and explicit build rules as needed)

bounce.exe: $(BOUNCE_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(BOUNCE_OBJECTS) $(BOUNCE_LIBS)

kraken.exe: $(KRAKEN_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(KRAKEN_OBJECTS) $(KRAKEN_LIBS)

krakenc.exe: $(KRAKENC_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(KRAKENC_OBJECTS) $(KRAKENC_LIBS)

InverseIterationMod.o inverseiterationmod.mod: InverseIterationMod.f90
	$(FC) -c $(FFLAGS) $<

RootFinderSecant.o rootfindersecant.mod: RootFinderSecantMod.f90
	$(FC) -c $(FFLAGS) $< 

KrakenMod.o krakenmod.mod: KrakenMod.f90
	$(FC) -c $(FFLAGS) $< 

KrakencMod.o krakencmod.mod: KrakencMod.f90
	$(FC) -c $(FFLAGS) $< 

BCImpedance.o: krakenmod.mod BCImpedance.f90
BCImpedancec.o: krakencmod.mod BCImpedancec.f90

bounce.o: krakencmod.mod bounce.f90
covar.o: covar.f90

kraken.o:  krakenmod.mod  kraken.f90
krakenc.o: krakencmod.mod krakenc.f90

kuping.o: kuping.f90
zbrentx.o: zbrentx.f90


