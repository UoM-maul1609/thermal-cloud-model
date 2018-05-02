BAM_DIR = bam

.PHONY: bam_code cleanall
CLEANDIRS = $(BAM_DIR) ./


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


wmicro_lib.a	:   nrtype.$(OBJ) nr.$(OBJ) nrutil.$(OBJ) locate.$(OBJ) polint.$(OBJ) \
				rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ) zbrent.$(OBJ) dfridr.$(OBJ) \
				hygfx.$(OBJ) microphysics.$(OBJ) advection_1d.$(OBJ) bam_code
	$(AR) rc wmicro_lib.a nrutil.$(OBJ) locate.$(OBJ) polint.$(OBJ) \
				rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ) zbrent.$(OBJ) dfridr.$(OBJ) \
				hygfx.$(OBJ) microphysics.$(OBJ) advection_1d.$(OBJ)
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
dfridr.$(OBJ)	: dfridr.f90
	$(FOR) dfridr.f90 $(FFLAGS)dfridr.$(OBJ)	
odeint.$(OBJ)	: odeint.f90
	$(FOR) odeint.f90 $(FFLAGS)odeint.$(OBJ)	
zbrent.$(OBJ)	: zbrent.f90
	$(FOR) zbrent.f90 $(FFLAGS2)zbrent.$(OBJ)	
advection_1d.$(OBJ) : advection_1d.f90 
	$(FOR) advection_1d.f90 $(FFLAGS)advection_1d.$(OBJ)
microphysics.$(OBJ) : microphysics.f90 advection_1d.$(OBJ) dfridr.$(OBJ) bam_code
	$(FOR) microphysics.f90 $(FFLAGS)microphysics.$(OBJ) -I$(BAM_DIR)
hygfx.$(OBJ) : hygfx.for 
	$(FOR) hygfx.for $(FFLAGS)hygfx.$(OBJ) 

bam_code:
	$(MAKE) -C $(BAM_DIR)

clean: 
	rm *.exe  *.o *.mod *~ \
	wmicro_lib.a

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done
	
