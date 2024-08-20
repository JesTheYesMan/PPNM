CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

decomp.dll: decomp.cs ../matrix/matrix.dll; $(MKLIB)

QRGS.dll: QRGS.cs linsol.cs ../matrix/matrix.dll; $(MKLIB) 
