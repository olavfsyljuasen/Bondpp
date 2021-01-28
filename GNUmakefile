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

# The envirmoment variable MYMAKEFILE contains the name of the appropriate Makefile.local
%    :  force
	cp $(MYMAKEFILE) Makefile.local 
	make -f makefile $*



force : ;

square_x64 : square.x64 squareforced.x64 squaresymmrot.x64 squaresymminv.x64 squaresymmaxis.x64 squaresymmdiag.x64

square_abel : square.abel squareforced.abel squaresymmrot.abel squaresymminv.abel squaresymmaxis.abel squaresymmdiag.abel

square_saga : square.saga squareforced.saga squaresymmrot.saga squaresymminv.saga squaresymmaxis.saga squaresymmdiag.saga


triangular_x64 : triangular.x64 triangularforced.x64 triangularsymmrot.x64 triangularsymminv.x64 triangularsymmaxis.x64 triangularsymmdiag.x64

triangular_abel : triangular.abel triangularforced.abel triangularsymmrot.abel triangularsymminv.abel triangularsymmaxis.abel triangularsymmdiag.abel

triangular_saga : triangular.saga triangularforced.saga triangularsymmrot.saga triangularsymminv.saga triangularsymmaxis.saga triangularsymmdiag.saga








