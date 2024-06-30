CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

default: Out.graph.svg Out.gnuplot.svg

Out1.gnuplot.svg: Out1.txt erf.data.txt Makefile
	echo '\
		set terminal $(subst .,,$(suffix $@)) background "white"; \
		set output "$@"; \
		set grid;\
		set key left;\
		set xlabel "x";\
		set ylabel "erf(x)";\
		set title "Error function";\
		plot [-4:4]\
			 "$<" with lines notitle \
			,"erf.data.txt" with points pointtype 6 notitle \
			,"erf.data.txt" using (-($$1)):(-($$2)) with points pointtype 6 not; \
	'| tee log.gpi | gnuplot

Out1.txt: main1.exe
	mono $< > $@

Out2.gnuplot.svg: Out2.txt fac.data.txt Makefile
	echo '\
		set terminal $(subst .,,$(suffix $@)) background "white"; \
		set output "$@"; \
		set grid;\
		set key left;\
		set xlabel "x";\
		set ylabel "gamma(x)";\
		set title "Gamma function";\
		set yrange [-30:30];\
		plot [-5:5.5]\
			"$<" with lines notitle \
			,"fac.data.txt" with points pointtype 6 notitle \
			,"fac.data.txt" using 1:2 with points pointtype 6 not; \
	'| tee log.gpi | gnuplot

Out2.txt: main2.exe
	mono $< > $@


Out3.gnuplot.svg: Out3.txt ln.fac.data.txt Makefile
	echo '\
		set terminal $(subst .,,$(suffix $@)) background "white"; \
		set output "$@"; \
		set grid;\
		set key left;\
		set xlabel "x";\
		set ylabel "lngamma(x)";\
		set title "ln Gamma function";\
		plot [-1:8]\
			"$<" with lines notitle \
			,"ln.fac.data.txt" with points pointtype 6 notitle \
			,"ln.fac.data.txt" using 1:2 with points pointtype 6 not; \
	'| tee log.gpi | gnuplot

Out3.txt: main3.exe
	mono $< > $@


main1.exe: main1.cs sfuns.dll ; $(MKEXE)
main2.exe: main2.cs sfuns.dll ; $(MKEXE)
main3.exe: main3.cs sfuns.dll ; $(MKEXE)
sfuns.dll : sfuns.cs ; $(MKLIB)

clean:
	$(RM) *.exe *.dll [Oo]ut* log*

