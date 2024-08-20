# Define the source code and libraries
CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # libraries
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

# Targets
all: Out.Harmonic_test.svg OutA.svg OutB.txt OutB.svg OutC.svg

Out.Harmonic_test.svg: data_A_Harmonic.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "t";\
		set grid;\
		set size 1,1;\
		set title "Testing harmonic";\
		plot "data_A_Harmonic.txt" using 1:2 with l t "Harmonic oscillator"\
		' | tee logA.gpi | gnuplot

# Create OutA.svg 
OutA.svg:  data_A_Oscilator.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "t";\
		set grid;\
		set size 1,1;\
		set title "Solution to Harmonic oscillator and damped pendulum";\
		plot "data_A_Oscilator.txt" using 1:2 with l t "theta(t)"\
		, "data_A_Oscilator.txt" using 1:3 with l t "omega(t)"\
		' | tee logA.gpi | gnuplot

# Create OutB.svg 
OutB.svg: data_B_Circular.txt data_B_Eliptical.txt data_B_GR.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "x";\
		set ylabel "y";\
		set size ratio -1;\
		set grid;\
		set xrange [-1.5:1.5];\
		set yrange [-1.5:2.5];\
		set term svg size 800, 400 ;\
		set multiplot layout 1,3 columns;\
		set title "Circular orbit";\
		plot "data_B_Circular.txt" using (1/$$2)*cos($$1):(1/$$2)*sin($$1) with lines notitle;\
		set title "Elliptic orbit";\
		plot "data_B_Eliptical.txt" using (1/$$2)*cos($$1):(1/$$2)*sin($$1) with lines notitle;\
		set title "Relativistic precession";\
		plot "data_B_GR.txt" using (1/$$2)*cos($$1):(1/$$2)*sin($$1) with lines notitle;\
		unset multiplot;\
		' | tee logB.gpi | gnuplot

# Create OutC.svg 
OutC.svg: data_C.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "x";\
		set ylabel "y";\
		set title "Figure-8 motion";\
		plot[-1.3:2] "$<" using 1:2 with lines t "Body A"\
		, "$<" using 3:4 with lines t "Body B"\
		, "$<" using 5:6 with lines t "Body C"\
		' | tee logC.gpi | gnuplot

# Generate data files
data_A_Harmonic.txt: mainA.exe
	mono $< 1> $@ 2> data_A_Oscilator.txt

data_B_Circular.txt: mainB.exe Makefile
	mono $< -u0:1 1> $@ 2>> OutB.txt

data_B_Eliptical.txt: mainB.exe Makefile
	mono $< -u0:1 -uprime0:-0.5 1> $@ 2>> OutB.txt

data_B_GR.txt: mainB.exe Makefile
	mono $< -u0:1 -uprime0:-0.5 -epsilon:0.01 -phiL:20 1> $@ 2>> OutB.txt

OutB.txt: 
	echo 'Test of ODE solver with no path-recording:' > $@

data_C.txt: mainC.exe
	mono $< > $@

# Build executables from C# source files
mainA.exe: mainA.cs genlib.dll ode.dll
	$(MKEXE)

mainB.exe: mainB.cs genlib.dll ode.dll
	$(MKEXE)

mainC.exe: mainC.cs genlib.dll ode.dll
	$(MKEXE)


genlib.dll: ../lib/matrix/vector.cs ../lib/genlist/genlist.cs \
	../lib/spline/spline.cs ../lib/matrix/matrix.cs \
	../lib/LinEq/QRGS.cs ../lib/LinEq/linsol.cs
	$(MKLIB)

ode.dll: ode.cs genlib.dll
	$(MKLIB)


clean:
	$(RM) *.dll *.exe [Oo]ut* [Dd]ata* [Ll]og*


