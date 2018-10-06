
EXECUTABLES = field.exe field3d.exe

FIELD_OBJECTS = field.o evaluate.o EvaluateADMod.o EvaluateCMMod.o ReadModes.o ../misc/libmisc.a
FIELD_LIBS = -L../misc -lmisc

FIELD3D_OBJECTS = field3d.o ElementMod.o Evaluate3DMod.o EvaluateGBMod.o evaluatepdq.o ReadModes.o ../misc/libmisc.a
FIELD3D_LIBS = -L../misc -lmisc

all: $(EXECUTABLES)
	@echo "KRAKEN field routines built"
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

field.exe: $(FIELD_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(FIELD_OBJECTS) $(FIELD_LIBS)

EvaluateCMMod.o evaluatecmmod.mod: EvaluateCMMod.f90
	$(FC) -c $(FFLAGS) $<

EvaluateADMod.o evaluateadmod.mod: EvaluateADMod.f90
	$(FC) -c $(FFLAGS) $<

field3d.exe: $(FIELD3D_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(FIELD3D_OBJECTS) $(FIELD3D_LIBS)

Evaluate3DMod.o evaluate3dmod.mod: Evaluate3DMod.f90
	$(FC) -c $(FFLAGS) $<

EvaluateGBMod.o evaluategbmod.mod: EvaluateGBMod.f90
	$(FC) -c $(FFLAGS) $<

ElementMod.o elementmod.mod: ElementMod.f90
	$(FC) -c $(FFLAGS) $<

ReadModes.o readmodes.mod: ReadModes.f90
	$(FC) -c $(FFLAGS) $<

evaluate.o: evaluate.f90
evaluatepdq.o: elementmod.mod evaluatepdq.f90

field.o:   readmodes.mod                evaluateadmod.mod evaluatecmmod.mod field.f90
field3d.o: readmodes.mod elementmod.mod evaluate3dmod.mod evaluategbmod.mod field3d.f90

