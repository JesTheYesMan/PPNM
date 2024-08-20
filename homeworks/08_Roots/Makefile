CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

all: OutA.txt OutBWF.svg OutBconv.svg OutC.txt

OutBconv.svg: dataBconv.txt Makefile
	echo '\
		set terminal svg background "white";\
		set out "$@";\
		set grid;\
		set ylabel "E";\
		f(x) = -0.5;\
		set term svg size 1200, 1200;\
		set multiplot layout 2,2 columns;\
		set title "r_{min} convergence";\
		set xlabel "r_{min}";\
		set xrange [0.01:0.01/2**9];\
		set logscale x 2;\
		set for [i=0:4] xtics (0.01,0.01/4**i);\
		plot "$<" using 1:2 with lp pt 7 notitle\
		,f(x) with lines lt rgb "green" title "exact value";\
		set title "r_{max} convergence";\
                set xlabel "r_{max}";\
		unset logscale x;\
		set for [i=0:9] xtics (6,6+i) ;\
		plot[6:16] "$<" using 3:4 with lp pt 7 notitle\
                ,f(x) with lines lt rgb "green" title "exact value";\
		set title "acc convergence";\
                set xlabel "acc";\
                set xrange [0.5:0.5/2**9];\
                set logscale x 2;\
                set for [i=0:4] xtics (0.5,0.5/4**i);\
                plot "$<" using 5:6 with lp pt 7 notitle\
                ,f(x) with lines lt rgb "green" title "exact value";\
		set title "eps convergence";\
                set xlabel "eps";\
                set xrange [0.5:0.5/2**9];\
                set logscale x 2;\
                set for [i=0:4] xtics (0.5,0.5/4**i);\
                plot "$<" using 7:8 with lp pt 7 notitle\
                ,f(x) with lines lt rgb "green" title "exact value";\
		unset multiplot;\
		'|tee logconv.gpi | gnuplot

OutBWF.svg: dataBWF.txt Makefile
	echo '\
                set terminal svg background "white" ;\
                set out "$@" ;\
                set xlabel "r";\
                set ylabel "Psi";\
                set xrange [0:8] ;\
                set title "n = 0 s-Wavefunction";\
                f1(x) = x*exp(-x);\
		plot "$<" using 1:2 with lines lt rgb "blue" title "Numerical Wavefunction" \
		,f1(x) with lines lt rgb "red"  title "Analytic Wavefunction"\
                '|tee logWF.gpi | gnuplot

OutA.txt: mainA.exe
	mono $< > $@

dataBWF.txt: mainB.exe
	mono $< 1> dataBconv.txt 2> $@

OutC.txt: mainC.exe
	mono $< > $@

mainA.exe: mainA.cs root.dll extlib.dll; $(MKEXE)

mainB.exe: mainB.cs root.dll extlib.dll; $(MKEXE) 

mainC.exe: mainC.cs root.dll extlib.dll; $(MKEXE)

extlib.dll: ../lib/matrix/vector.cs ../lib/matrix/matrix.cs ../lib/LinEq/QRGS.cs \
	../lib/ODE/ode.cs ../lib/genlist/genlist.cs; $(MKLIB)

root.dll: root.cs extlib.dll; $(MKLIB)

clean:
	$(RM) *.dll *.exe [Oo]ut* [Ll]og*
