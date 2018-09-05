FCC = gfortran

DM = ../MODULE/
DL = ../LIBRARY/

# Include files
IF1 = ../INCLUDE/para_geometry
IF2 = ../INCLUDE/para.h

a.out : g_module.o f_module.o fpara_module.o main.o cg.o clover.o \
	propa.o su3.o rand.o boundary_fermi.o initconf.o \
	measN.o \
	boundary_gauge.o lib_fermion.o
	$(FCC) g_module.o f_module.o fpara_module.o main.o cg.o clover.o \
	propa.o su3.o rand.o boundary_fermi.o initconf.o \
	measN.o \
	boundary_gauge.o lib_fermion.o

main.o : main.f f_module.o fpara_module.o $(IF1) $(IF2)
	$(FCC) -c main.f

measN.o : measN.f
	$(FCC) -c measN.f

initconf.o : $(DL)initconf.f
	$(FCC) -c $(DL)initconf.f

propa.o : propa.f f_module.o $(IF1)
	$(FCC) -c propa.f

cg.o : CGwithClover.f f_module.o fpara_module.o $(IF1) 
	$(FCC) -c -o cg.o CGwithClover.f

clover.o : clover.f f_module.o fpara_module.o $(IF1) 
	$(FCC) -c clover.f

su3.o : $(DL)su3.f g_module.o $(IF1) 
	$(FCC) -c $(DL)su3.f

boundary_gauge.o : $(DL)boundary_gauge.f g_module.o  $(IF1) $(IF2) 
	$(FCC) -c $(DL)boundary_gauge.f

boundary_fermi.o : $(DL)boundary_fermi.f f_module.o fpara_module.o \
                   $(IF1) $(IF2) 
	$(FCC) -c $(DL)boundary_fermi.f

lib_fermion.o : $(DL)lib_fermion.f f_module.o $(IF1)
	$(FCC) -c $(DL)lib_fermion.f

fpara_module.o : $(DM)fpara_module.f
	$(FCC) -c $(DM)fpara_module.f

f_module.o : $(DM)f_module.f g_module.o fpara_module.o $(IF1) $(IF2)
	$(FCC) -c $(DM)f_module.f

g_module.o : $(DM)g_module0.f $(IF1) 
	$(FCC) -c $(DM)g_module0.f -o g_module.o

rand.o : $(DL)rand.f
	$(FCC) -c $(DL)rand.f

clean :
	/bin/rm -f *.o *.mod $(DL)*.mod *.log
