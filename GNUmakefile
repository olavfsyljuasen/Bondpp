DATADIR= /scratch/sylju
HOMEDIR= /mn/felt/u1/sylju
ABELDIR= $(HOME)
SAGADIR= $(HOME)
OSXDIR = $(HOME)

%.x64 : force
	cp Makefile.x64 Makefile.local
	@echo `git describe`
	make -f makefile $*.exec

%.osx : force
	cp Makefile.osx Makefile.local
	@echo `git describe`
	make -f makefile $*.exec


%.abel : force
	cp Makefile.abel Makefile.local
	@echo `git describe`
	make -f makefile $*.exec

%.saga : force
	cp Makefile.saga Makefile.local
	@echo `git describe`
	make -f makefile $*.exec

%.fox : force
	cp Makefile.fox Makefile.local
	@echo `git describe`
	make -f makefile $*.exec

%.nbi : force
	cp Makefile.nbi Makefile.local
	@echo `git describe`
	make -f makefile $*.exec

# The envirmoment variable MYMAKEFILE contains the name of the appropriate Makefile.local
%    :  force
	cp $(MYMAKEFILE) Makefile.local 
	make -f makefile $*



force : ;

phonons_x64 : squareHeisenberg.x64 squareHeisenbergnoelastic.x64 squareHeisenbergelastic.x64 squareHeisenbergelasticclamped0.x64 squareHeisenbergelasticclamped1.x64 squareHeisenbergelasticclamped2.x64 

phonons_fox : squareHeisenberg.fox squareHeisenbergnoelastic.fox squareHeisenbergelastic.fox squareHeisenbergelasticclamped0.fox squareHeisenbergelasticclamped1.fox squareHeisenbergelasticclamped2.fox 

phonons_saga : squareHeisenberg.saga squareHeisenbergnoelastic.saga squareHeisenbergelastic.saga squareHeisenbergelasticclamped0.saga squareHeisenbergelasticclamped1.saga squareHeisenbergelasticclamped2.saga 



CrI3_x64: CrI3iso.x64 CrI3ssym.x64 CrI3nooffd.x64 CrI3.x64

CrI3_saga: CrI3iso.saga CrI3ssym.saga CrI3nooffd.saga CrI3.saga

#square_x64 : square.x64 squareiso.x64 squarenooffd.x64
#square_saga : square.saga squareiso.saga squarenooffd.saga

#squareXY_x64 : squareXY.x64 squareXYiso.x64 squareXYnooffd.x64
#squareXY_saga : squareXY.saga squareXYiso.saga squareXYnooffd.saga


square_x64:	squareHeisenberg.x64 squareHeisenbergelasticonly.x64 squareHeisenbergphononsonly.x64 squareHeisenbergphononsonly_print.x64 squareHeisenberg_print.x64 squareHeisenbergX.x64 squareHeisenbergelasticonlyX.x64 squareHeisenbergphononsonlyX.x64 

square_fox:	squareHeisenberg.fox squareHeisenbergelasticonly.fox squareHeisenbergphononsonly.fox squareHeisenbergphononsonly_print.fox squareHeisenberg_print.fox squareHeisenbergX.fox squareHeisenbergelasticonlyX.fox squareHeisenbergphononsonlyX.fox 

square_saga:	squareHeisenberg.saga squareHeisenbergelasticonly.saga squareHeisenbergphononsonly.saga squareHeisenbergphononsonly_print.saga squareHeisenberg_print.saga squareHeisenbergX.saga squareHeisenbergelasticonlyX.saga squareHeisenbergphononsonlyX.saga 


triangular_x64:	triangularHeisenberg.x64 triangularHeisenbergclamped.x64 triangularHeisenbergelasticonly.x64 triangularHeisenbergphononsonly.x64 triangularHeisenbergphononsonly_print.x64 triangularHeisenberg_print.x64 triangularHeisenbergX.x64 triangularHeisenbergelasticonlyX.x64 triangularHeisenbergphononsonlyX.x64 triangularHeisenbergX_print.x64

triangular_fox:	triangularHeisenberg.fox triangularHeisenbergclamped.fox triangularHeisenbergelasticonly.fox triangularHeisenbergphononsonly.fox triangularHeisenbergphononsonly_print.fox triangularHeisenberg_print.fox triangularHeisenbergX.fox triangularHeisenbergelasticonlyX.fox triangularHeisenbergphononsonlyX.fox triangularHeisenbergX_print.fox

triangular_saga: triangularHeisenberg.saga triangularHeisenbergclamped.saga triangularHeisenbergelasticonly.saga triangularHeisenbergphononsonly.saga triangularHeisenbergphononsonly_print.saga triangularHeisenberg_print.saga triangularHeisenbergX.saga triangularHeisenbergelasticonlyX.saga triangularHeisenbergphononsonlyX.saga triangularHeisenbergX_print.saga


triangularelasticclamped_x64:	triangularHeisenbergelastic0only.x64 triangularHeisenbergelastic1only.x64  triangularHeisenbergelastic2only.x64 triangularHeisenbergelastic01only.x64 triangularHeisenbergelastic02only.x64 triangularHeisenbergelastic12only.x64

triangularelasticclamped_fox:	triangularHeisenbergelastic0only.fox triangularHeisenbergelastic1only.fox  triangularHeisenbergelastic2only.fox triangularHeisenbergelastic01only.fox triangularHeisenbergelastic02only.fox triangularHeisenbergelastic12only.fox

triangularelasticclamped_saga:	triangularHeisenbergelastic0only.saga triangularHeisenbergelastic1only.saga  triangularHeisenbergelastic2only.saga triangularHeisenbergelastic01only.saga triangularHeisenbergelastic02only.saga triangularHeisenbergelastic12only.saga


squareelasticclamped_x64:	squareHeisenbergelastic0only.x64 squareHeisenbergelastic1only.x64  squareHeisenbergelastic2only.x64 squareHeisenbergelastic01only.x64 squareHeisenbergelastic02only.x64 squareHeisenbergelastic12only.x64

squareelasticclamped_fox:	squareHeisenbergelastic0only.fox squareHeisenbergelastic1only.fox  squareHeisenbergelastic2only.fox squareHeisenbergelastic01only.fox squareHeisenbergelastic02only.fox squareHeisenbergelastic12only.fox

squareelasticclamped_saga:	squareHeisenbergelastic0only.saga squareHeisenbergelastic1only.saga  squareHeisenbergelastic2only.saga squareHeisenbergelastic01only.saga squareHeisenbergelastic02only.saga squareHeisenbergelastic12only.saga


squareall_x64: squareHeisenbergP.x64 squareHeisenbergPE0.x64 squareHeisenbergPE1.x64 squareHeisenbergPE2.x64 squareHeisenbergPE01.x64 squareHeisenbergPE02.x64 squareHeisenbergPE12.x64 squareHeisenbergPE012.x64 squareHeisenbergE0.x64 squareHeisenbergE1.x64 squareHeisenbergE2.x64 squareHeisenbergE01.x64 squareHeisenbergE02.x64 squareHeisenbergE12.x64 squareHeisenbergE012.x64 squareHeisenbergE012noselfenergy.x64 squareHeisenbergE1noselfenergy.x64 

squareall_fox: squareHeisenbergP.fox squareHeisenbergPE0.fox squareHeisenbergPE1.fox squareHeisenbergPE2.fox squareHeisenbergPE01.fox squareHeisenbergPE02.fox squareHeisenbergPE12.fox squareHeisenbergPE012.fox squareHeisenbergE0.fox squareHeisenbergE1.fox squareHeisenbergE2.fox squareHeisenbergE01.fox squareHeisenbergE02.fox squareHeisenbergE12.fox squareHeisenbergE012.fox squareHeisenbergE012noselfenergy.fox squareHeisenbergE1noselfenergy.fox 

squareall_saga: squareHeisenbergP.saga squareHeisenbergPE0.saga squareHeisenbergPE1.saga squareHeisenbergPE2.saga squareHeisenbergPE01.saga squareHeisenbergPE02.saga squareHeisenbergPE12.saga squareHeisenbergPE012.saga squareHeisenbergE0.saga squareHeisenbergE1.saga squareHeisenbergE2.saga squareHeisenbergE01.saga squareHeisenbergE02.saga squareHeisenbergE12.saga squareHeisenbergE012.saga squareHeisenbergE012noselfenergy.saga squareHeisenbergE1noselfenergy.saga 

squareall_nbi: squareHeisenbergP.nbi squareHeisenbergPE0.nbi squareHeisenbergPE1.nbi squareHeisenbergPE2.nbi squareHeisenbergPE01.nbi squareHeisenbergPE02.nbi squareHeisenbergPE12.nbi squareHeisenbergPE012.nbi squareHeisenbergE0.nbi squareHeisenbergE1.nbi squareHeisenbergE2.nbi squareHeisenbergE01.nbi squareHeisenbergE02.nbi squareHeisenbergE12.nbi squareHeisenbergE012.nbi squareHeisenbergE012noselfenergy.nbi squareHeisenbergE1noselfenergy.nbi 


triangularall_x64: triangularHeisenbergP.x64 triangularHeisenbergPE0.x64 triangularHeisenbergPE1.x64 triangularHeisenbergPE2.x64 triangularHeisenbergPE01.x64 triangularHeisenbergPE02.x64 triangularHeisenbergPE12.x64 triangularHeisenbergPE012.x64 triangularHeisenbergE0.x64 triangularHeisenbergE1.x64 triangularHeisenbergE2.x64 triangularHeisenbergE01.x64 triangularHeisenbergE02.x64 triangularHeisenbergE12.x64 triangularHeisenbergE012.x64

triangularall_fox: triangularHeisenbergP.fox triangularHeisenbergPE0.fox triangularHeisenbergPE1.fox triangularHeisenbergPE2.fox triangularHeisenbergPE01.fox triangularHeisenbergPE02.fox triangularHeisenbergPE12.fox triangularHeisenbergPE012.fox triangularHeisenbergE0.fox triangularHeisenbergE1.fox triangularHeisenbergE2.fox triangularHeisenbergE01.fox triangularHeisenbergE02.fox triangularHeisenbergE12.fox triangularHeisenbergE012.fox

triangularall_saga: triangularHeisenbergP.saga triangularHeisenbergPE0.saga triangularHeisenbergPE1.saga triangularHeisenbergPE2.saga triangularHeisenbergPE01.saga triangularHeisenbergPE02.saga triangularHeisenbergPE12.saga triangularHeisenbergPE012.saga triangularHeisenbergE0.saga triangularHeisenbergE1.saga triangularHeisenbergE2.saga triangularHeisenbergE01.saga triangularHeisenbergE02.saga triangularHeisenbergE12.saga triangularHeisenbergE012.saga

triangularall_nbi: triangularHeisenbergP.nbi triangularHeisenbergPE0.nbi triangularHeisenbergPE1.nbi triangularHeisenbergPE2.nbi triangularHeisenbergPE01.nbi triangularHeisenbergPE02.nbi triangularHeisenbergPE12.nbi triangularHeisenbergPE012.nbi triangularHeisenbergE0.nbi triangularHeisenbergE1.nbi triangularHeisenbergE2.nbi triangularHeisenbergE01.nbi triangularHeisenbergE02.nbi triangularHeisenbergE12.nbi triangularHeisenbergE012.nbi















