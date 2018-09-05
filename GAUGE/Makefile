FCC = gfortran 
PFLAG =
OpFLAG =

DM = ../MODULE/
DL = ../LIBRARY/

MG = ../MODULE/g_module2.f
# g_module1.f : Wilson   Action
# g_module2.f : Improved Action

# Include files
IF1 = ../INCLUDE/para_geometry
IF2 = ../INCLUDE/para.h

a.out : g_module.o aniso_module.o main.o update.o meas1.o meas2.o \
	su3.o initconf.o saveconf.o rand.o boundary_gauge.o 
	$(FCC) $(PFLAG) $(OpFLAG)  g_module.o aniso_module.o main.o update.o meas1.o meas2.o \
	su3.o initconf.o saveconf.o rand.o boundary_gauge.o

main.o : main.f g_module.o aniso_module.o $(IF1) $(IF2)
	$(FCC)  $(PFLAG) $(OpFLAG) -c main.f

g_module.o : $(MG) $(IF1)
	$(FCC) $(PFLAG) $(OpFLAG)  -c $(MG) -o g_module.o

aniso_module.o : $(DM)aniso_module.f
	$(FCC) $(PFLAG) $(OpFLAG)  -c $(DM)aniso_module.f

update.o : update.f g_module.o aniso_module.o $(IF1) 
	$(FCC) $(PFLAG) $(OpFLAG)  -c update.f

meas1.o : meas1.f g_module.o $(IF1) 
	$(FCC) $(PFLAG) $(OpFLAG)  -c meas1.f

meas2.o : meas2.f g_module.o $(IF1) 
	$(FCC) $(PFLAG) $(OpFLAG)  -c meas2.f

su3.o : $(DL)su3.f g_module.o $(IF1) 
	$(FCC) $(PFLAG) $(OpFLAG)  -c $(DL)su3.f

initconf.o : $(DL)initconf.f
	$(FCC) $(PFLAG) $(OpFLAG)  -c $(DL)initconf.f

saveconf.o : $(DL)saveconf.f
	$(FCC) $(PFLAG) $(OpFLAG)  -c $(DL)saveconf.f

rand.o : ../LIBRARY/rand.f
	$(FCC) $(PFLAG) $(OpFLAG)  -c $(DL)rand.f

boundary_gauge.o : $(DL)boundary_gauge.f g_module.o $(IF1) $(IF2) 
	$(FCC) $(OpFLAG)  -c $(DL)boundary_gauge.f

# In case of debugging
# debug_g.o : $(DL)debug_g.f g_module.o $(IF1)
#	$(FCC) -c $(DL)debug_g.f

clean :
	/bin/rm -f *.o *.mod a.out
