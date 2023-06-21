OSNF_DIR = osnf

.PHONY: osnf_code cleanall
CLEANDIRS = $(OSNF_DIR) ./


DEBUG = -fbounds-check -g
MPI    =#-DMPI1
OPT    =-O3


FOR = gfortran -c  
FOR2 = gfortran  

AR = ar 
RANLIB = ranlib 
OBJ = o
FFLAGS = $(OPT) -w $(DEBUG)  -o 
FFLAGS2 =  $(DEBUG) -w -O3 -o 
VAR_TYPE = 1 # 0 single, 1 double


micro_lib.a	:   microphysics.$(OBJ) advection_1d.$(OBJ)
	$(AR) rc micro_lib.a \
				microphysics.$(OBJ) advection_1d.$(OBJ) \
				$(OSNF_DIR)/numerics.$(OBJ) $(OSNF_DIR)/zeroin.$(OBJ) $(OSNF_DIR)/sfmin.$(OBJ) \
                $(OSNF_DIR)/fmin.$(OBJ) $(OSNF_DIR)/r1mach.$(OBJ) \
                $(OSNF_DIR)/d1mach.$(OBJ) $(OSNF_DIR)/dfsid1.$(OBJ) \
                $(OSNF_DIR)/poly_int.$(OBJ) $(OSNF_DIR)/find_pos.$(OBJ) \
                $(OSNF_DIR)/svode.$(OBJ) \
                $(OSNF_DIR)/slinpk.$(OBJ) $(OSNF_DIR)/vode.$(OBJ) \
                $(OSNF_DIR)/dlinpk.$(OBJ) $(OSNF_DIR)/vode_integrate.$(OBJ) \
                $(OSNF_DIR)/erfinv.$(OBJ) $(OSNF_DIR)/tridiagonal.$(OBJ) \
                $(OSNF_DIR)/hygfx.$(OBJ) $(OSNF_DIR)/random.$(OBJ)				
advection_1d.$(OBJ) : advection_1d.f90 osnf_code
	$(FOR) advection_1d.f90 $(FFLAGS)advection_1d.$(OBJ) -I$(OSNF_DIR)
microphysics.$(OBJ) : microphysics.f90 advection_1d.$(OBJ) osnf_code
	$(FOR) microphysics.f90 $(FFLAGS)microphysics.$(OBJ) -I$(OSNF_DIR)

osnf_code:
	$(MAKE) -C $(OSNF_DIR)

clean: 
	rm *.exe  *.o *.mod *~ \
	micro_lib.a

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done

