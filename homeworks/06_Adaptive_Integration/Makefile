CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

all: OutA.txt OutA.Erf.svg OutB.txt OutC.txt

OutA.Erf.svg: dataA.txt dataAERF.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set key right bottom ;\
		set xlabel "x";\
		set ylabel "erf(x)";\
		set title "Error function";\
		set grid ;\
		plot [-3:3]\
			"$<" using 1:2 with lines t "Plots exercise approx" \
			,"$<" using 1:3 with lines t "Numerical integration" \
			,"dataAERF.txt" using 1:2 with points pointtype 6 notitle\
			,"dataAERF.txt" using (-($$1)):(-($$2)) with points pointtype 6 notitle ;\
	'| tee log.gpi | gnuplot 

OutA.txt: mainA.exe
	mono $< 1> $@ 2> dataA.txt

OutB.txt: mainB.exe
	mono $< > $@

OutC.txt: mainC.exe
	mono $< > $@

mainA.exe: mainA.cs int.dll extlib.dll; $(MKEXE)

mainB.exe: mainB.cs int.dll extlib.dll; $(MKEXE)

mainC.exe: mainC.cs int.dll extlib.dll; $(MKEXE)

extlib.dll: ../lib/matrix/vector.cs; $(MKLIB)

int.dll: int.cs; $(MKLIB)

clean:
	$(RM) *.dll *.exe [Oo]ut* [Ll]og*
