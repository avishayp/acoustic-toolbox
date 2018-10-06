
EXECUTABLES = bellhop.exe bellhop3d.exe

BELLHOP_OBJECTS   = bellhop.o   angleMod.o ArrMod.o bdryMod.o   bellhopMod.o ReadEnvironment.o sspMod.o influence.o ../misc/libmisc.a
BELLHOP_LIBS = -L../misc -lmisc

BELLHOP3D_OBJECTS = bellhop3D.o angleMod.o ArrMod.o bdry3DMod.o bellhopMod.o ReadEnvironment.o sspMod.o influence.o influence3D.o ../misc/libmisc.a
BELLHOP3D_LIBS = -L../misc -lmisc

all: $(EXECUTABLES)
	@echo "Bellhop built"
	@echo "*************"
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

bellhop.exe: $(BELLHOP_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(BELLHOP_OBJECTS) $(BELLHOP_LIBS)

bellhop3d.exe: $(BELLHOP3D_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(BELLHOP3D_OBJECTS) $(BELLHOP3D_LIBS)

angleMod.o anglemod.mod: angleMod.f90
	$(FC) -c $(FFLAGS) $< 

ArrMod.o arrmod.mod: ArrMod.f90
	$(FC) -c $(FFLAGS) $< 

bdryMod.o bdrymod.mod: bdryMod.f90
	$(FC) -c $(FFLAGS) $< 

bdry3DMod.o bdry3dmod.mod: bdry3DMod.f90
	$(FC) -c $(FFLAGS) $< 

bellhopMod.o bellhopmod.mod: bellhopMod.f90
	$(FC) -c $(FFLAGS) $<  

sspMod.o sspmod.mod: sspMod.f90
	$(FC) -c $(FFLAGS) $< 

bellhop.o:   arrmod.mod bellhopmod.mod anglemod.mod bdrymod.mod   sspmod.mod bellhop.f90
bellhop3D.o: arrmod.mod bellhopmod.mod anglemod.mod bdry3dmod.mod sspmod.mod bellhop3D.f90

ReadEnvironment.o:    anglemod.mod sspmod.mod ReadEnvironment.f90

influence.o:   influence.f90
influence3D.o: influence3D.f90


