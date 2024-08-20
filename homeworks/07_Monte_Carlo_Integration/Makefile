CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

all: OutA.txt OutA.svg OutB.txt OutB.svg OutC.txt

OutA.svg: dataA.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "N";\
		set ylabel "Errors";\
		set grid ;\
		f(x) = b + a*x**c;\
		b = 1; a = 90; c = -0.5;\
		set term svg size 1200 1200;\
		set multiplot layout 3,1 columns;\
		set title "Area of 2d circle";\
		fit f(x) "$<" using 1:3 via a,b;\
		plot "$<" using 1:2 with points t "Est. error"\
		,"$<" using 1:3 with points t "Exact error"\
		,f(x) with l title sprintf("fit: %.f + %.2f*N^{%.3f}",b,a,c);\
                set title "Volume of 3d half ball";\
                fit f(x) "$<" using 1:5 via a,b;\
                plot "$<" using 1:4 with points t "Est. error"\
                ,"$<" using 1:5 with points t "Exact error"\
                ,f(x) with l title sprintf("fit: %.f + %.2f*N^{%.3f}",b,a,c);\
                set title "Measure of 4d half ball";\
                fit f(x) "$<" using 1:7 via a,b;\
                plot "$<" using 1:6 with points t "Est. error"\
                ,"$<" using 1:7 with points t "Exact error"\
                ,f(x) with l title sprintf("fit: %.f + %.2f*N^{%.3f}",b,a,c);\
		unset multiplot;\
		'| tee logA.gpi | gnuplot

OutB.svg: dataB.txt dataA.txt Makefile
	echo '\
                set terminal svg background "white" ;\
                set out "$@" ;\
                set xlabel "N";\
                set ylabel "Errors";\
                set grid ;\
                set term svg size 1200 1200;\
                set multiplot layout 3,1 columns;\
                set title "Area of 2d circle";\
                plot "$<" using 1:2 with l t "Est. error"\
                ,"$<" using 1:3 with l t "Exact error"\
                ,"dataA.txt" using 1:3 with l t "Error of plainMC";\
                set title "Volume of 3d half ball";\
                plot "$<" using 1:4 with l t "Est. error"\
                ,"$<" using 1:5 with l t "Exact error"\
                ,"dataA.txt" using 1:5 with l t "Error of plainMC";\
                set title "Measure of 4d half ball";\
                plot "$<" using 1:6 with l t "Est. error"\
                ,"$<" using 1:7 with l t "Exact error"\
		,"dataA.txt" using 1:7 with l t "Error of plainMC";\
		unset multiplot;\
                '| tee logB.gpi | gnuplot

OutA.txt: mainA.exe
	mono $< 1> $@ 2> dataA.txt

OutB.txt: mainB.exe
	mono $< 1> $@ 2> dataB.txt

OutC.txt: mainC.exe
	mono $< > $@

mainA.exe: mainA.cs int.dll extlib.dll; $(MKEXE)

mainB.exe: mainB.cs int.dll extlib.dll; $(MKEXE)

mainC.exe: mainC.cs int.dll extlib.dll; $(MKEXE)

extlib.dll: ../lib/matrix/vector.cs ../lib/matrix/matrix.cs; $(MKLIB)

int.dll: int.cs extlib.dll; $(MKLIB)

clean: 
	$(RM) *.dll *.exe [Oo]ut* [Dd]ata* [Ll]og* *.log
