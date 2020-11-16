
PLANE_STRESS_SRC = src/precision.f90 \
	src/quicksort.f90 \
	src/plane_stress.f90

F90_FLAGS=
F90=gfortran

default: ${PLANE_STRESS_SRC}
	f2py -L/usr/lib -llapack -lblas ${F90_FLAGS} -c -m analysis_lib ${PLANE_STRESS_SRC}

debug: F90_FLAGS=--noopt --debug --f90flags="-fcheck=all"
debug: default

clean:
	rm -f *.o *.so
