CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)
TIME = time --output=$@ --append --format "$$N %e %U"

default: Out.txt Out.times.svg Makefile

Out.times.svg: outtimes.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "Matrix size N";\
		set ylabel "Decomposition time";\
		set title "Running time as function of size of matrices";\
		f(x)=b+(x/a)**c;\
		b=1; a=80; c=3;\
		fit f(x) "$<" via a,b;\
		plot "$<" with lp \
		,f(x) title sprintf("fit: (n/%.f)^{%.2f}+%.3f",a,c,b) \
		'|tee log.gpi |gnuplot

outtimes.txt: time.exe Makefile
	>$@
	for N in $$(seq 100 20 600); do \
		$(TIME) mono $< -msize:$$N 1>out 2>err ;\
	done ; wait

Out.txt: main.exe
	mono $< > $@

time.exe: time.cs matrix.dll QRGS.dll; $(MKEXE)

main.exe: main.cs matrix.dll QRGS.dll; $(MKEXE)

matrix.dll: ../lib/matrix/matrix.cs ../lib/matrix/vector.cs; $(MKLIB)

QRGS.dll: QRGS.cs matrix.dll; $(MKLIB)

clean: 
	$(RM) *.dll *.exe [Oo]ut* [Ll]og* *.[Ll]og
