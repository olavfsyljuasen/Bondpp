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















