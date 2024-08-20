CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

spline.dll: spline.cs ../matrix/matrix.cs ../matrix/vector.cs \
	../QRDecomp/QRGS.cs ../genlist/genlist.cs ../QRDecomp/linsol.cs; $(MKLIB)
