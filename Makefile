#
# To install the Acoustics Toolbox:
#
# 1) Uncomment the appropriate lines below to select your FORTRAN compiler
#    (also be sure to comment out all of the lines corresponding to the other compilers).

# 2) If you're using gfortran check the -march switch that selects the chip you're using.
#    Usually -march=native works

# 3) From a command line shell, run:
#    % make clean
#    % make

# on some machines you need to say -mcmodel=medium (or large) to allow for variables larger than 2 gig

# *** Windows ***
# If you don't have a FORTRAN compiler, the MinGW ("Minimalist GNU for Windows") compiler suite
# is a good choice as it is much easier to install than Cygwin. See http://www.mingw.org/ for details.
# Not sure if this is still necessary, but we used to have to change the options -O3 to -O2 below when using gfortran

# *** Linux ***
# Most Linux distributions have gfortran already packaged and it can be installed with the respective 
# package manager (e.g. apt-get, dnf, yum). The packaged versions of the LAPACK library are generally
# compatible with Krakel. If you want statically linked executables, "-static" works with gfortran.

# *** Mac ***
# The option to create an executable that uses a static library has been a problem.
# The gcc community seems to feel that dynamic vs. static libraries is an issue of religious importance.
# I prefer static libraries because users have a lot fewer problems installing the code when they don't have
# to worry about getting the LD_LIBRARY_PATH set, or the fact that some programs (e.g. Matlab) may reset
# that path to an incompatible version.
# -static is supposed to work if gfortran has been compiled to enable that option.
# However, crt0.o comes up as misssing
# -Bstatic worked

# ______________________________________________________________________________

# *** ifort
# These lines are used under Mac OSX; the syntax is different under Windows
# You also need to use xiar instead of ar
# choose the best architecture (target machine) using the -x switch

# need -heap-arrays to avoid stack overflows for big runs ...
export FC=ifort
# export FFLAGS= -O3 -parallel -axSSE4.2 -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I../misc
export FFLAGS= -fast              -parallel -axAVX -nologo              -inline-level=2 -assume byterecl -threads -heap-arrays -I../misc
export FFLAGS= -O3 -funroll-loops -parallel -no-prec-div -axAVX -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I../misc
export FFLAGS= -O3 -funroll-loops -parallel -xAVX -no-prec-div -axAVX -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I../misc
export FFLAGS= -O3 -xHost -qopt-report -funroll-loops -parallel -no-prec-div -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I../misc
export FFLAGS= -O3 -xHost -funroll-loops -no-prec-div -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I../misc

# export FFLAGS= -O3 -fast -ipo     -parallel -no-prec-div                                -assume byterecl -heap-arrays -I../misc   # used on Polyhedron.com site

# -parallel off for export:
# export FFLAGS= -O3 -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I../misc
# export FFLAGS= -O3 -ipo -funroll-loops -xAVX -no-prec-div -axAVX -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I../misc

# compilation diagnostics on:
export FFLAGS= -O2 -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I../misc -check -traceback

# runtime diagnostic on as well:
export FFLAGS= -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I../misc -check all -ftrapuv -fpe0 -gen-interfaces -traceback
#export FFLAGS= -nologo -O0 -inline-level=2 -assume byterecl -threads -heap-arrays -I../misc -check all -ftrapuv -gen-interfaces -traceback
#export FFLAGS= -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I../misc -check all -ftrapuv -traceback

# profiling:
# export FFLAGS= -g -O3 -xHost -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I../misc -traceback

# ______________________________________________________________________________

# *** GNU Compiler Collection GFORTRAN

# use -march=generic if you get warning messages about instructions that don't make sense
# -march=generic assumes an old Intel architecture that the newer versions can all execute (slowly)
# -march=native should normally be the best; however, it produced AVX instructions on the Mac that the default assembler could not process
# -O2 was the highest level of optimization that worked under Windows
#
# -static can be used to tell gfortran not to rely on a dynamic link library (the compiler may or may not support)
# -static does not seem to work on Macs though, and produces larger executables
# Have had various problems where some installed dynamic link library is incompatible with the one the compiler used and expects at run time
# For instance, Matlab changes paths and may point to an incompatible library.
# One user found that it was necessary to delete /usr/local/gfortran/lib/libquadmath.dylib to force a static link. See:
# http://stackoverflow.com/questions/17590525/correct-way-to-statically-link-in-gfortran-libraries-on-osx
#
# The -Wa,-q flag can be used to select the Mac CLANG assembler instead of the GNU assembler
# At one time that was necessary to get the AVX operations; however, I saw no speed benefit
# -march=corei7-avx works on my Mac

export FC=gfortran

# export FFLAGS= -march=native -Wall -std=gnu -O2 -I../misc
# export FFLAGS= -mtune=generic -Wall -std=gnu -O3 -I../misc
# export FFLAGS= -march=native -Wall -std=gnu -O3 -I../misc -ffast-math -funroll-all-loops -msse3 -fomit-frame-pointer -mtune=native -Q
# export FFLAGS= -march=corei7 -Bstatic -Waliasing -Wampersand -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -std=f2008 -O3 -ffast-math -funroll-all-loops -fomit-frame-pointer -mtune=native -I../misc
# -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk/usr/lib

#export FFLAGS= -march=corei7 -Bstatic -Waliasing -Wampersand -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -std=f2008 -O3 -ffast-math -funroll-all-loops -fomit-frame-pointer -I../misc

# export FFLAGS= -march=native -Waliasing -Wampersand -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wa,-q -std=f2008 -O3 -ffast-math -funroll-all-loops -fomit-frame-pointer -mtune=native -I../misc

export FFLAGS= -march=native -Bstatic -Waliasing -Wampersand -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -std=gnu -O3 -ffast-math -funroll-all-loops -fomit-frame-pointer -mtune=native -I../misc

#compilation and run-time diagnostics on:
#export FFLAGS= -march=corei7 -ffpe-trap=invalid,zero,overflow -Wall -std=gnu -O1 -fcheck=all -fbacktrace -I../misc
#export FFLAGS= -march=native -ffpe-trap=zero,overflow -Wall -std=gnu -O1 -fcheck=all -I../misc

# profiling:
# I read that the -pg flag is needed for profiling, but there's some problem with the library and it doesn't compile
# It does compile with just the -p and -g flags but does not appear to have enough info for the xcode instruments
#export FFLAGS= -p -g -pg -march=native -Wall -std=f2008 -I../misc
#export FFLAGS= -p -g -march=native -Bstatic -Wa,-q -Waliasing -Wampersand -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -std=f2008 -O3 -ffast-math -funroll-all-loops -fcheck=all -I../misc

# ______________________________________________________________________________

# *** g95
# this is no longer working with the Acoustics Toolbox
# export FC=g95
# export FFLAGS= -Wall -std=f2003 -O3 -I../misc

# compilation diagnostics on:
# export FFLAGS= -Wall -std=f2003 -ftrace=full -fbounds-check -I../misc
# export FFLAGS = -pg -std=f2003 -I../misc

# ______________________________________________________________________________

# *** Portland Group FORTRAN
# -Mnoframe caused erroneous results
# -Munroll  caused erroneous results
# these are defaults under -fast, so can't use -fast either
# export FC=pgfortran
# export FFLAGS= -Mconcur -I../misc
# export FFLAGS= -fast -I../misc
# export FFLAGS= -O2 -Munroll=c:1 -Mnoframe -Mlre -Mpre -Mvect=sse -Mcache_align -Mflushz -Mvect -I../misc
# export FFLAGS= -O2 -Mlre -Mpre -Mvect=sse -Mcache_align -Mflushz -Mvect -I../misc

# compilation diagnostics on:
# export FFLAGS= -g -Minfo=ccff -Minform=inform -C -I../misc

# ______________________________________________________________________________

export RM=rm
export CC=gcc
export CFLAGS=-g

# KRAKEL is commented out below because it requires the LAPACK library.
# If you have the LAPACK library installed on your system, first edit the
# LAPACK_LIBS variable below to point to your installation, then you can
# uncomment the make commands below (uncomment both the "all" and "install"
# targets).

export LAPACK_LIBS = -llapack

all:
	(cd misc;	make -k all)
	(cd tslib;	make -k all)
	(cd Bellhop;	make -k all)
	(cd Kraken;	make -k all)
	(cd KrakenField;	make -k all)
	#(cd Krakel;	make -k all)
	(cd Scooter;	make -k all)
	@echo " "
	@echo "***********************************"
	@echo "***** Acoustics Toolbox built *****"
	@echo "***********************************"

#install:
#	(cd Bellhop;	make -k install)
#	(cd Kraken;	make -k install)
#	(cd KrakenField;        make -k install)
	#(cd Krakel;	make -k install)
#	(cd Scooter;	make -k install)
#	@echo " "
#	@echo "***************************************"
#	@echo "***** Acoustics Toolbox installed *****"
#	@echo "***************************************"

clean:
	-rm -f bin/*.exe
	(cd misc;	make -k -i clean)
	(cd tslib;	make -k -i clean)
	(cd Bellhop;	make -k -i clean)
	(cd Kraken;	make -k -i clean)
	(cd KrakenField;	make -k -i clean)
	(cd Krakel;	make -k -i clean)
	(cd Scooter;	make -k -i clean)
	(cd tests;	make -k -i clean)

