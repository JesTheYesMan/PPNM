CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)
TIME = time --output=$@ --append --format "$$N %e %U"

all: Out_A.txt Out_B.txt main_A.dll

Out_A.txt: main_A.exe 
	mono $< > $@

Out_B.txt: main_B.exe 
	mono $< > $@


Ns := $(shell seq 100 25 300)

main_A.exe: main_A.cs mat.dll; $(MKEXE)

main_B.exe: main_B.cs main_A.dll mat.dll; $(MKEXE)

hydrogen.exe: hydrogen.cs mat.dll; $(MKEXE)

mat.dll: ../homeworks/lib/LinEq/QRGS.cs ../homeworks/lib/minimization/min.cs ../homeworks/lib/matrix/matrix.cs ../homeworks/lib/matrix/vector.cs jacobi.cs; $(MKLIB) 

main_A.dll: main_A.cs mat.dll; $(MKLIB)

clean: 
	$(RM) *.dll *.exe [Oo]ut* [Ll]og* *.[Ll]og
