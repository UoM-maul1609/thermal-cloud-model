MPM_DIR = mpm
WMM_DIR = wmm
BAM_DIR = wmm/bam
SFVT_DIR = sfvt

.PHONY: mpm_code cleanall
.PHONY: wmm_code cleanall
.PHONY: sfvt_code cleanall
CLEANDIRS = $(MPM_DIR) $(WMM_DIR) $(WMM_DIR)/bam $(SFVT_DIR) ./


DEBUG = -fbounds-check -g
MPI    =#-DMPI1
OPT    =-O3

# these three lines should be edited for your system. On systems 
# that do not have separate fortran and c libraries, set NETCDF_FOR and NETCDF_C
# to the same, and set NETCDF_LIB to -lnetcdf (i.e. without the extra f)
#NETCDF_FOR=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.4-mac/
#NETCDF_C=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.1.1-mac/
NETCDF_LIB=-lnetcdff 

NETCDFLIB=-L ${NETCDF_FOR}/lib/  \
          -L ${NETCDF_C}/lib/
NETCDFMOD= ${NETCDF_FOR}/include/


FOR = gfortran -c  
FOR2 = gfortran  

AR = ar 
RANLIB = ranlib 
OBJ = o
FFLAGS = $(OPT)  $(DEBUG)  -o 
FFLAGS2 =  $(DEBUG) -O3 -o 


main.exe	:  model_lib.a  main.$(OBJ) variables.$(OBJ) initialisation.$(OBJ) \
			 driver_code.$(OBJ) thermal_code.$(OBJ) io_code.$(OBJ) advection_2d.$(OBJ) \
			 advection_1d.$(OBJ) mpm_code wmm_code sfvt_code
	$(FOR2) $(FFLAGS2)main.exe main.$(OBJ) variables.$(OBJ) initialisation.$(OBJ) \
	     driver_code.$(OBJ) thermal_code.$(OBJ) io_code.$(OBJ) advection_2d.$(OBJ) \
	     		 $(MPM_DIR)/micro_lib.a $(WMM_DIR)/wmicro_lib.a $(BAM_DIR)/bam_lib.a \
	     		  $(SFVT_DIR)/model_lib.a \
	     		 -lm model_lib.a \
		${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)
model_lib.a	:   nrtype.$(OBJ) nr.$(OBJ) nrutil.$(OBJ) locate.$(OBJ) polint.$(OBJ) \
				rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ) zbrent.$(OBJ) \
				hygfx.$(OBJ) 
	$(AR) rc model_lib.a nrutil.$(OBJ) locate.$(OBJ) polint.$(OBJ) \
				rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ) zbrent.$(OBJ) \
				hygfx.$(OBJ) 
locate.$(OBJ)	: locate.f90
	$(FOR) locate.f90 $(FFLAGS)locate.$(OBJ)
polint.$(OBJ)	: polint.f90
	$(FOR) polint.f90 $(FFLAGS)polint.$(OBJ)
nrtype.$(OBJ)	: nrtype.f90
	$(FOR) nrtype.f90 $(FFLAGS)nrtype.$(OBJ)
nr.$(OBJ)	: nr.f90 
	$(FOR) nr.f90 $(FFLAGS)nr.$(OBJ)
nrutil.$(OBJ)	: nrutil.f90
	$(FOR) nrutil.f90 $(FFLAGS)nrutil.$(OBJ)
rkqs.$(OBJ)	: rkqs.f90
	$(FOR) rkqs.f90 $(FFLAGS)rkqs.$(OBJ)	
rkck.$(OBJ)	: rkck.f90
	$(FOR) rkck.f90 $(FFLAGS)rkck.$(OBJ)	
odeint.$(OBJ)	: odeint.f90
	$(FOR) odeint.f90 $(FFLAGS)odeint.$(OBJ)	
zbrent.$(OBJ)	: zbrent.f90
	$(FOR) zbrent.f90 $(FFLAGS2)zbrent.$(OBJ)	
variables.$(OBJ) : variables.f90 
	$(FOR) variables.f90 $(FFLAGS)variables.$(OBJ)
initialisation.$(OBJ) : initialisation.f90 mpm_code
	$(FOR) initialisation.f90 $(FFLAGS)initialisation.$(OBJ) -I$(MPM_DIR)
driver_code.$(OBJ) : driver_code.f90 thermal_code.$(OBJ) io_code.$(OBJ) \
		advection_2d.$(OBJ) advection_1d.$(OBJ) mpm_code wmm_code sfvt_code
	$(FOR) driver_code.f90 $(FFLAGS)driver_code.$(OBJ) -I$(MPM_DIR) -I$(WMM_DIR) -I$(SFVT_DIR)
thermal_code.$(OBJ) : thermal_code.f90 
	$(FOR) thermal_code.f90 $(FFLAGS)thermal_code.$(OBJ)
io_code.$(OBJ) : io_code.f90 
	$(FOR) io_code.f90 -I ${NETCDFMOD} $(FFLAGS)io_code.$(OBJ)
advection_1d.$(OBJ) : advection_1d.f90 
	$(FOR) advection_1d.f90 $(FFLAGS)advection_1d.$(OBJ)
advection_2d.$(OBJ) : advection_2d.f90 
	$(FOR) advection_2d.f90 $(FFLAGS)advection_2d.$(OBJ)
hygfx.$(OBJ) : hygfx.for 
	$(FOR) hygfx.for $(FFLAGS)hygfx.$(OBJ) 
main.$(OBJ)   : main.f90 variables.$(OBJ) initialisation.$(OBJ) driver_code.$(OBJ) \
			 thermal_code.$(OBJ) io_code.$(OBJ) advection_2d.$(OBJ) advection_1d.$(OBJ) \
			 mpm_code wmm_code
	$(FOR)  main.f90 -I ${NETCDFMOD} -I${BAM_DIR} $(FFLAGS)main.$(OBJ) 
	

mpm_code:
	$(MAKE) -C $(MPM_DIR)

wmm_code:
	$(MAKE) -C $(WMM_DIR)

sfvt_code:
	$(MAKE) -C $(SFVT_DIR)

clean: 
	rm *.exe  *.o *.mod *~ \
	model_lib.a

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done
	
	
