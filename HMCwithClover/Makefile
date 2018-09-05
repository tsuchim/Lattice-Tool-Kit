FCC = gfortran
PFLAG = 
OpFLAG =

DM = ../MODULE/
DL = ../LIBRARY/
DF = ../FERMIONwithClover/
DG = ../GAUGE/

# Include files
IF1 = ../INCLUDE/para_geometry
#IF2 = ../INCLUDE/para.h
IF2 = 

a.out  : fpara_module.o aniso_module.o g_module.o f_module.o hmc_module.o \
	main.o hmc_clover.o cg.o clover.o hmc.o staple.o meas1.o \
	boundary_gauge.o boundary_fermi.o su3.o initconf.o saveconf.o \
	su3_hmc.o meas.o lib_fermion.o rand.o debug1.o 
	$(FCC) $(PFLAG) $(OpFLAG)  fpara_module.o aniso_module.o \
	g_module.o f_module.o hmc_module.o \
	main.o hmc_clover.o cg.o clover.o hmc.o staple.o meas1.o \
	boundary_gauge.o  boundary_fermi.o su3.o initconf.o \
	saveconf.o su3_hmc.o meas.o lib_fermion.o rand.o debug1.o

main.o : main.f g_module.o f_module.o hmc_module.o fpara_module.o \
         aniso_module.o $(IF1) $(IF2) \
         debug1.o
	$(FCC) $(PFLAG) $(OpFLAG) -c main.f

hmc_clover.o : hmc_clover.f g_module.o f_module.o hmc_module.o \
	fpara_module.o $(IF1) debug1.o
	$(FCC) $(PFLAG) $(OpFLAG) -c hmc_clover.f

hmc.o : hmc.f g_module.o f_module.o hmc_module.o fpara_module.o $(IF1) \
	debug1.o
	$(FCC) $(PFLAG) $(OpFLAG) -c hmc.f

meas.o : meas.f g_module.o $(IF1) 
	$(FCC) $(PFLAG) $(OpFLAG) -c meas.f

meas1.o : $(DG)meas1.f g_module.o $(IF1) 
	$(FCC) $(PFLAG) $(OpFLAG) -c $(DG)meas1.f

staple.o : staple.f g_module.o aniso_module.o $(IF1) 
	$(FCC) $(PFLAG) $(OpFLAG) -c staple.f

su3_hmc.o : $(DL)su3_hmc.f g_module.o hmc_module.o $(IF1)
	$(FCC) $(PFLAG) $(OpFLAG) -c $(DL)su3_hmc.f

su3.o : $(DL)su3.f g_module.o $(IF1)
	$(FCC) $(PFLAG) $(OpFLAG) -c $(DL)su3.f

initconf.o : $(DL)initconf.f
	$(FCC) -c $(DL)initconf.f

saveconf.o : $(DL)saveconf.f
	$(FCC) -c $(DL)saveconf.f

lib_fermion.o : $(DL)lib_fermion.f f_module.o $(IF1) 
	$(FCC) $(PFLAG) $(OpFLAG) -c $(DL)lib_fermion.f

cg.o : $(DF)CGwithClover.f fpara_module.o f_module.o $(IF1) 
	$(FCC) $(PFLAG) $(OpFLAG) -c -o cg.o $(DF)CGwithClover.f

clover.o : $(DF)clover.f fpara_module.o f_module.o $(IF1) 
	$(FCC) $(PFLAG) $(OpFLAG) -c $(DF)clover.f

boundary_gauge.o : $(DL)boundary_gauge.f g_module.o $(IF1) $(IF2) 
	$(FCC) $(PFLAG) $(OpFLAG) -c $(DL)boundary_gauge.f

boundary_fermi.o : $(DL)boundary_fermi.f f_module.o fpara_module.o \
                   $(IF1) $(IF2) 
	$(FCC) $(PFLAG) $(OpFLAG) -c $(DL)boundary_fermi.f

fpara_module.o : $(DM)fpara_module.f
	$(FCC) $(PFLAG) $(OpFLAG) -c $(DM)fpara_module.f

aniso_module.o : $(DM)aniso_module.f
	$(FCC) -c $(DM)aniso_module.f

g_module.o : $(DM)g_module0.f $(IF1) 
	$(FCC) $(PFLAG) $(OpFLAG) -c $(DM)g_module0.f -o g_module.o

f_module.o : $(DM)f_module.f g_module.o fpara_module.o $(IF1) $(IF2) 
	$(FCC) $(PFLAG) $(OpFLAG) -c $(DM)f_module.f

hmc_module.o : $(DM)hmc_module.f g_module.o f_module.o $(IF1) $(IF2) 
	$(FCC) $(PFLAG) $(OpFLAG) -c $(DM)hmc_module.f

rand.o : $(DL)rand.f
	$(FCC) -c $(DL)rand.f

# In case of debugging
debug1.o : $(DM)debug1.f g_module.o hmc_module.o $(IF1)
	$(FCC) -c $(DM)debug1.f

clean :
	/bin/rm -f *.o *.mod $(DL)*mod $(DF)*mod $(DG)*mod
