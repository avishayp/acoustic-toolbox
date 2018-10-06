
EXECUTABLES = fields.exe scooter.exe sparc.exe

FIELDS_OBJECTS = fields.o TransformG.o pade.o ../misc/libmisc.a
FIELDS_LIBS = -L../tslib -lts -L../misc -lmisc

SCOOTER_OBJECTS = scooter.o factorNPZ.o ScooterMod.o ../misc/libmisc.a
SCOOTER_LIBS = -L../tslib -lts -L../misc -lmisc

SPARC_OBJECTS = sparc.o factorNP.o SparcMod.o ../misc/libmisc.a
SPARC_LIBS = -L../tslib -lts -L../misc -lmisc

all: $(EXECUTABLES)
	@echo "Scooter and Sparc built"
	@echo "***********************"
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

fields.exe: $(FIELDS_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(FIELDS_OBJECTS) $(FIELDS_LIBS)

scooter.exe: $(SCOOTER_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(SCOOTER_OBJECTS) $(SCOOTER_LIBS)

sparc.exe: $(SPARC_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(SPARC_OBJECTS) $(SPARC_LIBS)

ScooterMod.o scootermod.mod: ScooterMod.f90
	$(FC) -c $(FFLAGS) $< 

SparcMod.o sparcmod.mod: SparcMod.f90
	$(FC) -c $(FFLAGS) $< 

factorNP.o: factorNP.f90
factorNPZ.o: factorNPZ.f90

fields.o: fields.f90
TransformG.o: TransformG.f90
pade.o: pade.f90

scooter.o: scootermod.mod scooter.f90
sparc.o: sparcmod.mod sparc.f90

