MPM_DIR = mpm
WMM_DIR = wmm
BAM_DIR = wmm/bam
BAM_DIR2 = pamm/bam
SFVT_DIR = sfvt
PAMM_DIR = pamm
OSNF_DIR = sfvt/osnf

.PHONY: mpm_code cleanall
.PHONY: wmm_code cleanall
.PHONY: sfvt_code cleanall
.PHONY: pamm_code cleanall
CLEANDIRS = $(MPM_DIR) $(WMM_DIR) $(WMM_DIR)/bam $(SFVT_DIR) $(PAMM_DIR) $(PAMM_DIR)/sfvt \
       $(PAMM_DIR)/sfvt/osnf $(PAMM_DIR)/bam/osnf $(PAMM_DIR)/bam  \
       $(MPM_DIR)/osnf $(WMM_DIR)/bam/osnf $(SFVT_DIR)/osnf ./


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
FFLAGS = $(OPT)  $(DEBUG) -w -o 
FFLAGS2 =  $(DEBUG) -w -O3 -o 
VAR_TYPE = 1 # 0 single, 1 double


main.exe	:  model_lib.a  main.$(OBJ) variables.$(OBJ) initialisation.$(OBJ) \
			 driver_code.$(OBJ) thermal_code.$(OBJ) io_code.$(OBJ) advection_2d.$(OBJ) \
			 advection_1d.$(OBJ) mpm_code wmm_code sfvt_code pamm_code
	$(FOR2) $(FFLAGS2)main.exe main.$(OBJ) variables.$(OBJ) initialisation.$(OBJ) \
	     driver_code.$(OBJ) thermal_code.$(OBJ) io_code.$(OBJ) advection_2d.$(OBJ) \
	     		 $(PAMM_DIR)/pmicro_lib.a $(MPM_DIR)/micro_lib.a $(WMM_DIR)/wmicro_lib.a \
	     		 -lm model_lib.a -L$(PAMM_DIR) \
		${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)
model_lib.a	:   sfvt_code 
	$(AR) rc model_lib.a $(OSNF_DIR)/numerics.$(OBJ) $(OSNF_DIR)/zeroin.$(OBJ) $(OSNF_DIR)/sfmin.$(OBJ) \
				$(OSNF_DIR)/fmin.$(OBJ) $(OSNF_DIR)/r1mach.$(OBJ) \
                $(OSNF_DIR)/d1mach.$(OBJ) $(OSNF_DIR)/dfsid1.$(OBJ) \
                $(OSNF_DIR)/poly_int.$(OBJ) $(OSNF_DIR)/find_pos.$(OBJ) \
                $(OSNF_DIR)/svode.$(OBJ) \
                $(OSNF_DIR)/slinpk.$(OBJ) $(OSNF_DIR)/vode.$(OBJ) \
                $(OSNF_DIR)/dlinpk.$(OBJ) $(OSNF_DIR)/vode_integrate.$(OBJ) \
                $(OSNF_DIR)/erfinv.$(OBJ) $(OSNF_DIR)/tridiagonal.$(OBJ) \
                $(OSNF_DIR)/hygfx.$(OBJ) $(OSNF_DIR)/random.$(OBJ)	
variables.$(OBJ) : variables.f90 sfvt_code
	$(FOR) variables.f90 -I$(OSNF_DIR) $(FFLAGS)variables.$(OBJ)
initialisation.$(OBJ) : initialisation.f90 mpm_code sfvt_code
	$(FOR) initialisation.f90 $(FFLAGS)initialisation.$(OBJ) -I$(MPM_DIR) -I$(OSNF_DIR) 
driver_code.$(OBJ) : driver_code.f90 thermal_code.$(OBJ) io_code.$(OBJ) \
		advection_2d.$(OBJ) advection_1d.$(OBJ) mpm_code wmm_code sfvt_code pamm_code
	$(FOR) driver_code.f90 $(FFLAGS)driver_code.$(OBJ) -I$(MPM_DIR) -I$(WMM_DIR) -I$(SFVT_DIR) \
	         -I$(PAMM_DIR) -I$(OSNF_DIR) 
thermal_code.$(OBJ) : thermal_code.f90 sfvt_code
	$(FOR) thermal_code.f90 $(FFLAGS)thermal_code.$(OBJ) -I$(OSNF_DIR) 
io_code.$(OBJ) : io_code.f90 sfvt_code
	$(FOR) io_code.f90 -I ${NETCDFMOD} $(FFLAGS)io_code.$(OBJ) -I$(OSNF_DIR) 
advection_1d.$(OBJ) : advection_1d.f90 sfvt_code
	$(FOR) advection_1d.f90 $(FFLAGS)advection_1d.$(OBJ) -I$(OSNF_DIR) 
advection_2d.$(OBJ) : advection_2d.f90 sfvt_code
	$(FOR) advection_2d.f90 $(FFLAGS)advection_2d.$(OBJ) -I$(OSNF_DIR) 
main.$(OBJ)   : main.f90 variables.$(OBJ) initialisation.$(OBJ) driver_code.$(OBJ) \
			 thermal_code.$(OBJ) io_code.$(OBJ) advection_2d.$(OBJ) advection_1d.$(OBJ) \
			 mpm_code wmm_code pamm_code
	$(FOR)  main.f90 -I ${NETCDFMOD} -I${MPM_DIR} -I${WMM_DIR} -I${PAMM_DIR} \
	    $(FFLAGS)main.$(OBJ) 
	

mpm_code:
	$(MAKE) -C $(MPM_DIR)

wmm_code:
	$(MAKE) -C $(WMM_DIR)

sfvt_code:
	$(MAKE) -C $(SFVT_DIR)

pamm_code:
	$(MAKE) -C $(PAMM_DIR)

clean: 
	rm *.exe  *.o *.mod *~ \
	model_lib.a

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done
	
	
