CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros for filtering .cs and .dll files
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)
TIME = time --output=$@ --append --format "$$N %e %U"

all: Out.svg outdata.txt

# f(x) = (cos(x)+sin(x)/2)*exp(-x/2); damped solution

Out.svg: outdata.txt dataTable.txt Makefile
	echo '\
                set terminal svg background "white" ;\
                set out "$@" ;\
                set xlabel "x";\
                set ylabel "g";\
		set key outside;\
		g(x) = sin(4*x-1)*exp(-x*x);\
		set term svg size 1000, 1200 ;\
		set multiplot layout 3,1 columns;\
                set title "Neural Network for interpolation";\
                set grid;\
                plot[-1.1:1.1] g(x) with l lt 0 t "g(x)"\
		,"dataTable.txt" with p pt 6 ps 0.5 t "Training data"\
		,"$<" using 1:2 with l t "ANN 4 nodes" linecolor rgb "#4627E2";\
		set title "derivative:";\
		plot[-1.1:1.1]"$<" using 1:3 with l t "ANN 4 nodes" linecolor rgb "#4627E2";\
		set title "integral from x=-1";\
		plot[-1.1:1.1]"$<" using 1:4 with l t "ANN 4 nodes" linecolor rgb "#4627E2";\
		unset multiplot;\
                '|tee log.gpi | gnuplot


outdata.txt: main.exe
	mono $< -output:dataTable.txt > $@

main.exe: main.cs ann.dll extlib.dll
	$(MKEXE)

extlib.dll: ../lib/matrix/vector.cs ../lib/minimization/min.cs \
	../lib/matrix/matrix.cs ../lib/LinEq/QRGS.cs \
	../lib/Integration/int.cs
	$(MKLIB)

ann.dll: ann.cs extlib.dll
	$(MKLIB)

clean:
	$(RM) *.exe *.dll Out* log*

