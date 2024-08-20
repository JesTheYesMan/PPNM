CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

all: extlib.dll root.dll

extlib.dll: ../matrix/vector.cs ../matrix/matrix.cs ../QRDecomp/QRGS.cs \
	../ODE/ode.cs ../genlist/genlist.cs; $(MKLIB)

root.dll: root.cs extlib.dll; $(MKLIB)
	
clean:
	$(RM) *.dll
