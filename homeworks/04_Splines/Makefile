CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

all: Out.Linear_interp.svg Out.Quadratic_interp.svg Out.Cubic_interp.svg \
	out_linear_1.txt out_linear_2.txt \
	outquad.txt dataquad.txt \
	outcub1.txt outcub2.txt

Out.Linear_interp.svg: out_linear_1.txt out_linear_2.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "x";\
		set ylabel "y";\
		set size 1,1;\
		set key outside;\
                set term svg size 600, 800 ;\
                set multiplot layout 2,1 columns ;\
		set title "Linear interpolation";\
		plot "out_linear_2.txt" with p pt 7 ps 0.5 title "Data"\
		,"$<" using 1:2 with l notitle;\
		set title "Integral of Linear interpolation";\
		plot "$<" using 1:3 with l notitle;\
		unset multiplot;\
		'|tee log.gpi | gnuplot

Out.Quadratic_interp.svg: outquad.txt dataquad.txt Makefile
	echo '\
                set terminal svg background "white" ;\
                set out "$@" ;\
                set xlabel "x";\
                set ylabel "y";\
		set size 1,1;\
                set term svg size 600, 1200 ;\
		set multiplot layout 3,1 columns;\
                set title "Quadratic spline";\
		plot "dataquad.txt" using 1:2 with p pt 7 ps 0.5 title "Data"\
		,"dataquad.txt" using 1:3 with p pt 7 ps 0.5 notitle \
                ,"dataquad.txt" using 1:4 with p pt 7 ps 0.5 notitle \
		,"$<" using 1:2 with l notitle linecolor rgb "#4627E2"\
                ,"$<" using 1:3 with l notitle linecolor rgb "#27E23A"\
                ,"$<" using 1:4 with l notitle linecolor rgb "#F0860C";\
		set title "Quadratic spline derivatives";\
                plot "$<" using 1:5 with l notitle linecolor rgb "#4627E2"\
                ,"$<" using 1:6 with l notitle linecolor rgb "#27E23A"\
                ,"$<" using 1:7 with l notitle linecolor rgb "#F0860C";\
		set title "Quadratic spline integrals";\
		plot "$<" using 1:8 with l notitle linecolor rgb "#4627E2"\
                ,"$<" using 1:9 with l notitle linecolor rgb "#27E23A"\
		,"$<" using 1:10 with l notitle linecolor rgb "#F0860C";;\
		unset multiplot;\
                '|tee log.gpi | gnuplot

Out.Cubic_interp.svg: outcub1.txt outcub2.txt Makefile
	echo '\
                set terminal svg background "white" ;\
                set out "$@" ;\
                set xlabel "x";\
                set ylabel "y";\
		set size 1,1;\
		set key outside;\
		set term svg size 600, 1200 ;\
                set multiplot layout 3,1 columns ;\
                set title "Cubic spline";\
                plot "outcub2.txt" with p pt 7 ps 0.5 title "Data"\
		,"" smooth csplines t "gnuplot spline"\
                ,"$<" using 1:2 with l title "Homemade spline";\
                set title "Cubic spline derivative";\
                plot "$<" using 1:3 with l notitle;\
		set title "Cubic spline integral";\
		plot "$<" using 1:4 with l notitle;\
                unset multiplot;\
                '|tee log.gpi | gnuplot

out_linear_1.txt: main_linear.exe
	mono $< 1> $@ 2> out_linear_2.txt

outquad.txt: main_quadratic.exe
	mono $< -output:dataquad.txt 2> $@ 1> OutqTest.txt


outcub1.txt: main_cubic.exe
	mono $< 1> $@ 2> outcub2.txt

main_quadratic.exe: main_quadratic.cs spline.dll extlib.dll; $(MKEXE)

main_linear.exe: main_linear.cs spline.dll extlib.dll; $(MKEXE)

main_cubic.exe: main_cubic.cs spline.dll extlib.dll; $(MKEXE)

spline.dll: spline.cs extlib.dll; $(MKLIB)
	
extlib.dll: ../lib/matrix/matrix.cs ../lib/matrix/vector.cs\
       	../lib/LinEq/QRGS.cs ../lib/LinEq/linsol.cs \
	../lib/genlist/genlist.cs; $(MKLIB)

clean: 
	$(RM) *.dll *.exe [Oo]ut* [Ll]og*
