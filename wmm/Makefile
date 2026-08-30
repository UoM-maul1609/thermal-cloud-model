BAM_DIR = bam
OSNF_DIR = bam/osnf

.PHONY: bam_code osnf cleanall
CLEANDIRS = $(BAM_DIR) $(OSNF_DIR) ./


DEBUG = -fbounds-check -g
MPI    =#-DMPI1
OPT    =-O3


FOR = gfortran -c  
FOR2 = gfortran  

AR = ar 
RANLIB = ranlib 
OBJ = o
FFLAGS = $(OPT)  $(DEBUG)  -o 
FFLAGS2 =  $(DEBUG) -O0 -o 
VAR_TYPE = 1 # 0 single, 1 double


wmicro_lib.a	:   microphysics.$(OBJ) advection_1d.$(OBJ) bam_code 
	$(AR) rc wmicro_lib.a microphysics.$(OBJ) advection_1d.$(OBJ) \
				$(OSNF_DIR)/numerics.$(OBJ) $(BAM_DIR)/bulk_activation_module.$(OBJ) \
				$(OSNF_DIR)/random.$(OBJ) 	
advection_1d.$(OBJ) : advection_1d.f90 bam_code
	$(FOR) advection_1d.f90 $(FFLAGS)advection_1d.$(OBJ)  -I$(OSNF_DIR)
microphysics.$(OBJ) : microphysics.f90 advection_1d.$(OBJ) bam_code
	$(FOR) microphysics.f90 $(FFLAGS)microphysics.$(OBJ) -I$(BAM_DIR) -I$(OSNF_DIR)

bam_code:
	$(MAKE) -C $(BAM_DIR)

clean: 
	rm *.exe  *.o *.mod *~ \
	wmicro_lib.a

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done
	
