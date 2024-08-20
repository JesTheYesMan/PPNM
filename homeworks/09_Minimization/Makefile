CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

all: Out_A.txt Out_B_plot.svg Out_data_B.txt 

Out_B_plot.svg: higgs_data.txt Out_data_B.txt Makefile
	echo '\
                set terminal svg background "white" ;\
                set out "$@" ;\
                set xlabel "E [GeV]";\
                set ylabel "{/symbol s} (background subtracted)";\
		set title "Resonance fit";\
		set grid;\
                plot "$<" with errorbars t "Data"\
		, "Out_data_B.txt" with l t "Breit-Wigner fit"\
                '|tee logB.gpi | gnuplot

Out_A.txt: main_A.exe
	mono $< > $@

Out_data_B.txt: main_B.exe
	mono $< < higgs_data.txt 1> $@ 2> logBhiggs.txt


main_A.exe: main_A.cs extlib.dll min.dll; $(MKEXE)

main_B.exe: main_B.cs extlib.dll min.dll; $(MKEXE)


extlib.dll: ../lib/matrix/matrix.cs ../lib/matrix/vector.cs \
	../lib/genlist/genlist.cs ../lib/LinEq/QRGS.cs; $(MKLIB)

min.dll: min.cs extlib.dll; $(MKLIB)

clean:
	$(RM) *.dll *.exe [Oo]ut* [Ll]og*

