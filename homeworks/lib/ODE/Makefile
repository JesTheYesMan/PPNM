CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

ode.dll: ode.cs ../../libs/matrix/vector.cs ../../libs/genlist/genlist.cs; $(MKLIB)
