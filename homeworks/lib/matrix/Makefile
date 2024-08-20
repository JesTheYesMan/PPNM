CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

Out.txt: main.exe
	mono $< > $@

main.exe: main.cs vector.dll matrix.dll; $(MKEXE)

vec.dll: vec.cs; $(MKLIB)

vector.dll: vector.cs; $(MKLIB)

matrix.dll: matrix.cs vector.cs; $(MKLIB)

clean: 
	$(RM) *.dll *.exe [Oo]ut*

