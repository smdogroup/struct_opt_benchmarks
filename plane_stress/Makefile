
PLANE_STRESS_SRC = src/precision.f90 \
	src/quicksort.f90 \
	src/plane_stress.f90

SOLID_SRC = src/precision.f90 \
	src/quicksort.f90 \
	src/solid.f90

F90_FLAGS=
F90=gfortran

default: plane_stress solid

plane_stress: ${PLANE_STRESS_SRC}
	f2py -L/usr/lib -llapack -lblas ${F90_FLAGS} -c -m analysis_lib ${PLANE_STRESS_SRC}

solid: ${SOLID_SRC}
	f2py -L/usr/lib -llapack -lblas ${F90_FLAGS} -c -m solid_lib ${SOLID_SRC}

debug: F90_FLAGS=--noopt --debug --f90flags="-fcheck=all"
debug: default

clean:
	rm -f *.o *.so
