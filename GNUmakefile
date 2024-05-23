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

square_x64 : square.x64 squareiso.x64 squarenooffd.x64
square_saga : square.saga squareiso.saga squarenooffd.saga

squareXY_x64 : squareXY.x64 squareXYiso.x64 squareXYnooffd.x64
squareXY_saga : squareXY.saga squareXYiso.saga squareXYnooffd.saga










