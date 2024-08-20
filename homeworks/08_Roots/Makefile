CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

all: OutA.txt

OutA.txt: mainA.exe
	mono $< > $@

dataBWF.txt: mainB.exe
	mono $< 1> dataBconv.txt 2> $@


mainA.exe: mainA.cs root.dll extlib.dll; $(MKEXE)


extlib.dll: ../lib/matrix/vector.cs ../lib/matrix/matrix.cs ../lib/LinEq/QRGS.cs \
	../lib/ODE/ode.cs ../lib/genlist/genlist.cs; $(MKLIB)

root.dll: root.cs extlib.dll; $(MKLIB)

clean:
	$(RM) *.dll *.exe [Oo]ut* [Ll]og*
