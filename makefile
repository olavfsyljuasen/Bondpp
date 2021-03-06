include Makefile.local


PPSOURCES =  bondpp.C bondpp.h matrixroutines.h overload.h bravaislattices.h vecmat.h modeldef.h couplings.h rules.h fourierplans.h phonons.h inputparameters.h makefile Makefile.local globalheader.h RunParameter.h rnddef.h


SQUARESOURCES      =  $(PPSOURCES) 
TRIANGULARSOURCES  =  $(PPSOURCES) triangularrules.h
HONEYCOMBSOURCES   =  $(PPSOURCES) honeycombrules.h
KAGOMESOURCES      =  $(PPSOURCES) kagomerules.h

SIMPLECUBICSOURCES =  $(PPSOURCES) simplecubicrules.h
DIAMONDSOURCES     =  $(PPSOURCES) diamondrules.h
PYROCHLORESOURCES  =  $(PPSOURCES) pyrochlorerules.h

all : planar 3d

planar : square.exec triangular.exec kagome.exec honeycomb.exec 

3d: simplecubic.exec diamond.exec pyrochlore.exec

sq : square.exec squaresymmrot.exec squaresymminv.exec



%.run  : %.exec
	./$<

%.exec : %.o 
	$(CCC) $(LINKOPTS) -o $@ $^ $(LIBDIR) $(LIBS)
	mv $*.exec $(BINDIR)/$*.`git describe`$(EXTENSION)



#square.o  : $(SQUARESOURCES)
#	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUARELATTICE  -D#QSPACEDIRECTLATTICE -DNSPINCOMPONENTS=3 -DMAKEHERMITIAN -DFFTS_INPLACE -c -o $@ $<

square.o  : $(SQUARESOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUAREPHONONS -DSQUARELATTICE  -DQSPACEDIRECTLATTICE -DPHONONS -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFFTS_INPLACE -c -o $@ $<

CrI3.o  : $(SQUARESOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DCrI3 -DSQUARELATTICE  -DQSPACEDIRECTLATTICE -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFFTS_INPLACE -c -o $@ $<




spotless:       
	make clean
	rm -f *.ps

clean   :
	rm -f *.dvi
	rm -f *.aux
	rm -f *.log
	rm -f *~
	rm -f core
	rm -f *.o	
	rm -f *.exec
	rm -f *.d

