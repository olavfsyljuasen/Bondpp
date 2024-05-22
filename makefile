include Makefile.local


SOURCES =  bondpp.C bondpp.h matrixroutines.h overload.h bravaislattices.h vecmat.h modeldef.h couplings.h modelcouplings.h rules.h phonons.h inputparameters.h symmetryroutines.h observables.h makefile Makefile.local globalheader.h RunParameter.h rnddef.h mynumbertypes.h



all : planar 3d

planar : square.exec triangular.exec kagome.exec honeycomb.exec 

3d: simplecubic.exec diamond.exec pyrochlore.exec

sq : square.exec squaresymmrot.exec squaresymminv.exec



%.run  : %.exec
	./$<

%.exec : %.o 
	$(CCC) $(LINKOPTS) -o $@ $^ $(LIBDIR) $(LIBS)
	mv $*.exec $(BINDIR)/$*.`git describe`$(EXTENSION)



#square.o  : $(SOURCES)
#	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUARELATTICE  -D#QSPACEDIRECTLATTICE -DNSPINCOMPONENTS=3 -DMAKEHERMITIAN -DFFTS_INPLACE -c -o $@ $<

#square.o  : $(SOURCES)
#	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUAREPHONONS -DSQUARELATTICE  -DQSPACEDIRECTLATTICE -DPHONONS -DXYDISPLACEMENTS -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY  -DFFTS_INPLACE -c -o $@ $<

square.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUARENOPHONONS -DSPINISOTROPIC -DSQUARELATTICE  -DQSPACEDIRECTLATTICE -DUSELASTSIGMA -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DNBRRANGE=2 -DFFTS_INPLACE -c -o $@ $<


squareHeisenberg.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUAREPHONONS -DSQUARELATTICE  -DQSPACEDIRECTLATTICE -DPHONONS -DXYDISPLACEMENTS -DUSELASTSIGMA -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DNBRRANGE=2 -DFFTS_INPLACE -c -o $@ $<

squareHeisenbergomit.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUAREPHONONS -DSQUARELATTICE  -DQSPACEDIRECTLATTICE -DPHONONS -DXYDISPLACEMENTS -DUSELASTSIGMA -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DNBRRANGE=2 -DOMITPHONONCONT -DFFTS_INPLACE -c -o $@ $<

squareHeisenberg_print.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUAREPHONONS -DSQUARELATTICE  -DQSPACEDIRECTLATTICE -DPHONONS -DXYDISPLACEMENTS -DUSELASTSIGMA -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DNBRRANGE=2 -DFFTS_INPLACE -DPRINTPHONONS -c -o $@ $<

squareHeisenberg_rnd.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUAREPHONONS -DSQUARELATTICE  -DQSPACEDIRECTLATTICE -DPHONONS -DXYDISPLACEMENTS -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DNBRRANGE=2 -DFFTS_INPLACE -c -o $@ $<

squareHeisenbergelastic.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUAREPHONONS -DSQUARELATTICE  -DQSPACEDIRECTLATTICE -DPHONONS -DXYDISPLACEMENTS -DUSELASTSIGMA -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DNBRRANGE=2 -DFFTS_INPLACE -DELASTICONLY -DOMITPHONONCONT -c -o $@ $<


squareHeisenbergelasticclamped0.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUAREPHONONS -DSQUARELATTICE  -DQSPACEDIRECTLATTICE -DPHONONS -DXYDISPLACEMENTS -DUSELASTSIGMA -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DNBRRANGE=2 -DFFTS_INPLACE -DELASTICONLY -DOMITPHONONCONT -DONEEPSILONCOMPONENTCLAMPED -DEPSILONCOMPONENTCLAMPED=0 -c -o $@ $<

squareHeisenbergelasticclamped1.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUAREPHONONS -DSQUARELATTICE  -DQSPACEDIRECTLATTICE -DPHONONS -DXYDISPLACEMENTS -DUSELASTSIGMA -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DNBRRANGE=2 -DFFTS_INPLACE -DELASTICONLY -DOMITPHONONCONT -DONEEPSILONCOMPONENTCLAMPED -DEPSILONCOMPONENTCLAMPED=1 -c -o $@ $<

squareHeisenbergelasticclamped2.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUAREPHONONS -DSQUARELATTICE  -DQSPACEDIRECTLATTICE -DPHONONS -DXYDISPLACEMENTS -DUSELASTSIGMA -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DNBRRANGE=2 -DFFTS_INPLACE -DELASTICONLY -DOMITPHONONCONT -DONEEPSILONCOMPONENTCLAMPED -DEPSILONCOMPONENTCLAMPED=2 -c -o $@ $<

squareHeisenbergnoelastic.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUAREPHONONS -DSQUARELATTICE  -DQSPACEDIRECTLATTICE -DPHONONS -DXYDISPLACEMENTS -DUSELASTSIGMA -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DNBRRANGE=2 -DFFTS_INPLACE -DNOELASTIC -c -o $@ $<


squareHeisenbergnophonons.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUAREPHONONS -DSQUARELATTICE  -DQSPACEDIRECTLATTICE -DPHONONS -DXYDISPLACEMENTS -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DNBRRANGE=2 -DFFTS_INPLACE -DOMITPHONONCONT -DOMITELASTICCONT -c -o $@ $<


squarelayeredHeisenberg.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUARELAYEREDPHONONS -DSQUARELATTICE  -DQSPACEDIRECTLATTICE -DPHONONS -DXYZDISPLACEMENTS -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DFFTS_INPLACE -c -o $@ $<

#cubicHeisenberg.o  : $(SOURCES)
#	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DCUBICPHONONS -DCUBICLATTICE -DQSPACEDIRECTLATTICE -DPHONONS -DXYZDISPLACEMENT#S -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY  -DFAKEHEISENBERG -DCPOSITIVE -DFFTS_INPLACE -c -o $@ $<

cubicHeisenberg.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DCUBICPHONONS -DCUBICLATTICE -DQSPACEDIRECTLATTICE -DPHONONS -DXYZDISPLACEMENTS -DUSELASTSIGMA -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DFFTS_INPLACE -c -o $@ $<

cubicHeisenbergnoselfenergy.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DCUBICPHONONS -DCUBICLATTICE -DQSPACEDIRECTLATTICE -DPHONONS -DXYZDISPLACEMENTS -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DNOSELFENERGY -DFAKEHEISENBERG -DCPOSITIVE -DFFTS_INPLACE -c -o $@ $<

cubicHeisenbergelasticonly.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DCUBICPHONONS -DCUBICLATTICE -DQSPACEDIRECTLATTICE -DPHONONS -DELASTICONLY -DXYZDISPLACEMENTS -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DFFTS_INPLACE -c -o $@ $<

cubicHeisenbergelastic.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DCUBICPHONONS -DCUBICLATTICE -DQSPACEDIRECTLATTICE -DPHONONS -DXYZDISPLACEMENTS -DUSELASTSIGMA -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DFFTS_INPLACE -DELASTICONLY -DOMITPHONONCONT -c -o $@ $<

cubicHeisenbergnoelastic.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DCUBICPHONONS -DCUBICLATTICE -DQSPACEDIRECTLATTICE -DPHONONS -DXYZDISPLACEMENTS -DUSELASTSIGMA -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DFFTS_INPLACE -DNOELASTIC -c -o $@ $<


cubicHeisenbergnophonons.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DCUBICLATTICE -DQSPACEDIRECTLATTICE -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY  -DFAKEHEISENBERG -DCPOSITIVE -DFFTS_INPLACE -c -o $@ $<

triangularHeisenberg.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DHEXAGONALBRAVAISLATTICE -DNSUBLATTICES=1 -DTRIANGULARPHONONS -DHEXAGONALLATTICE  -DQSPACEDIRECTLATTICE -DPHONONS -DXYDISPLACEMENTS -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFAKEHEISENBERG -DCPOSITIVE -DFFTS_INPLACE -c -o $@ $<

triangularlayeredHeisenberg.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DHEXAGONALBRAVAISLATTICE -DNSUBLATTICES=1 -DTRIANGULARLAYEREDPHONONS -DHEXAGONALLATTICE -DQSPACEDIRECTLATTICE -DPHONONS -DXYZDISPLACEMENTS -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY  -DFAKEHEISENBERG -DCPOSITIVE -DFFTS_INPLACE -c -o $@ $<

fccHeisenberg.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DFCCBRAVAISLATTICE -DNSUBLATTICES=1 -DFCCPHONONS -DFCCLATTICE -DQSPACEDIRECTLATTICE -DPHONONS -DXYZDISPLACEMENTS -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY  -DFAKEHEISENBERG -DCPOSITIVE -DFFTS_INPLACE -c -o $@ $<

bccHeisenberg.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DBCCBRAVAISLATTICE -DNSUBLATTICES=1 -DBCCPHONONS -DBCCLATTICE -DQSPACEDIRECTLATTICE -DPHONONS -DXYZDISPLACEMENTS -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY  -DFAKEHEISENBERG -DCPOSITIVE -DFFTS_INPLACE -c -o $@ $<



CrI3nooffd.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DHEXAGONALBRAVAISLATTICE -DNSUBLATTICES=2 -DCrI3 -DHEXAGONALLATTICE -DQSPACEDIRECTLATTICE -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFFTS_INPLACE -DFFTS_INPLACE -DNORANDOMSPINOFFDIAGONALS -c -o $@ $<

CrI3ssym.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DHEXAGONALBRAVAISLATTICE -DNSUBLATTICES=2 -DCrI3 -DHEXAGONALLATTICE -DQSPACEDIRECTLATTICE -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFFTS_INPLACE -DFFTS_INPLACE -DNORANDOMSPINOFFDIAGONALS -DEQUALRANDOMSPINDIAGONALS -c -o $@ $<


CrI3.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DHEXAGONALBRAVAISLATTICE -DNSUBLATTICES=2 -DCrI3 -DHEXAGONALLATTICE -DQSPACEDIRECTLATTICE -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFFTS_INPLACE -DFFTS_INPLACE -c -o $@ $<


CrI3iso.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DHEXAGONALBRAVAISLATTICE -DNSUBLATTICES=2 -DCrI3 -DHEXAGONALLATTICE -DQSPACEDIRECTLATTICE -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFFTS_INPLACE -DFFTS_INPLACE -DSPINISOTROPIC -c -o $@ $<

squareiso.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUARENOPHONONS -DSQUARELATTICE -DQSPACEDIRECTLATTICE -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFFTS_INPLACE -DFFTS_INPLACE -DSPINISOTROPIC -c -o $@ $<

#square.o  : $(SOURCES)
#	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUARENOPHONONS -DSQUARELATTICE -DQSPACEDIRECTLATTICE -DRANDOMINITIALIZATION -DFFTS_INPLACE -c -o $@ $<

squarenooffd.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUARENOPHONONS -DSQUARELATTICE -DQSPACEDIRECTLATTICE -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFFTS_INPLACE -DNORANDOMSPINOFFDIAGONALS -c -o $@ $<

squareXY.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUAREXY -DSQUARELATTICE -DQSPACEDIRECTLATTICE -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFFTS_INPLACE -DFFTS_INPLACE -c -o $@ $<

squareXYiso.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUAREXY -DSQUARELATTICE -DQSPACEDIRECTLATTICE -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFFTS_INPLACE -DFFTS_INPLACE -DSPINISOTROPIC -c -o $@ $<

squareXYnooffd.o  : $(SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBICBRAVAISLATTICE -DNSUBLATTICES=1 -DSQUAREXY -DSQUARELATTICE -DQSPACEDIRECTLATTICE -DRANDOMINITIALIZATION -DFORCEINVERSIONSYMMETRY -DFFTS_INPLACE -DFFTS_INPLACE -DNORANDOMSPINOFFDIAGONALS -c -o $@ $<


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

