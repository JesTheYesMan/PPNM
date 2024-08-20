CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

all: extlib.dll min.dll

extlib.dll: ../matrix/matrix.cs ../matrix/vector.cs \
	../genlist/genlist.cs ../LinEq/QRGS.cs; $(MKLIB)

min.dll: min.cs extlib.dll; $(MKLIB)
