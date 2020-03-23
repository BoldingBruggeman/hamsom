#SASCII start run from physics climatology
#SATLAS start run biology from WOCA ATlas values
#else remember defining the input file in the code
#.SUFFIXES: .f .F .F90 .f90 .o .mod
#COMP=mpiifort -g -cpp -O2 -fp-model source -assume byterecl -DMPI -DMPIP
COMP=mpiifort -g -C -traceback -check bounds -cpp -fp-model source -assume byterecl -mcmodel=medium -DMPIP -DMPI
OTHERFLAGS = $(shell /sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-intel14/bin/nf-config --fflags)
LIBFLAGS = $(shell /sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-intel14/bin/nf-config --flibs)
OBJ=ncio.o hydrodynamics_fabm.o  ice_dynamics.o  setup_fabm.o  subs.o  trmice_cd2swr.o linear.o parallel_fabm.o main_fabm.o
ECOSMO_ttt.out :$(OBJ)
	$(COMP) $(OPT) $(OBJ) $(LIBFLAGS) -o ECOSMO_ttt.out
clean: 
	rm $(OBJ) 
ncio.o: ncio.f90
	$(COMP) $(OTHERFLAGS) -c ./ncio.f90
hydrodynamics_fabm.o: hydrodynamics_fabm.F
	$(COMP) $(OPT) -c hydrodynamics_fabm.F

ice_dynamics.o: ice_dynamics.F
	$(COMP) $(OPT) -c ice_dynamics.F

main_fabm.o: main_fabm.F
	$(COMP) $(OPT) $(OTHERFLAGS) -c  main_fabm.F

setup_fabm.o : setup_fabm.F
	$(COMP) $(OPT) -c setup_fabm.F

subs.o : subs.F
	$(COMP) $(OPT) -c subs.F

trmice_cd2swr.o: trmice_cd2swr.F
	$(COMP) $(OPT) -c trmice_cd2swr.F

parallel_fabm.o: parallel_fabm.F
	$(COMP) $(OPT) -c parallel_fabm.F

linear.o: linear.F
	$(COMP) $(OPT) -c linear.F

co2_dyn.o: co2_dyn.F90
	$(COMP) $(OPT) -c co2_dyn.F90


