CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

all: Out.Linear_interp.svg Out.Quadratic_interp.svg Out.Cubic_interp.svg \
	out_A_data_linear_1.txt out_A_data_linear_2.txt \
	out_B_data_quad.txt dataquad.txt \
	out_C_data_cubic_1.txt out_C_data_cubic_2.txt

Out.Linear_interp.svg: out_A_data_linear_1.txt out_A_data_linear_2.txt Makefile
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
		plot "out_A_data_linear_2.txt" with p pt 7 ps 0.5 title "Data"\
		,"$<" using 1:2 with l notitle;\
		set title "Integral of Linear interpolation";\
		plot "$<" using 1:3 with l notitle;\
		unset multiplot;\
		'|tee log.gpi | gnuplot

Out.Quadratic_interp.svg: out_B_data_quad.txt dataquad.txt Makefile
	echo '\
	set terminal svg background "white";\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set size 1,1;\
	set term svg size 600,800;\
	set multiplot layout 3,1 columns;\
	set title "Quadratic spline";\
	plot "dataquad.txt" using 1:2 with points pointtype 7 pointsize 0.5 title "Data",\
	"dataquad.txt" using 1:3 with points pointtype 7 pointsize 0.5 notitle,\
	"dataquad.txt" using 1:4 with points pointtype 7 pointsize 0.5 notitle,\
	"$<" using 1:2 with lines notitle linecolor rgb "#008080",\
	"$<" using 1:3 with lines notitle linecolor rgb "#DC143C",\
	"$<" using 1:4 with lines notitle linecolor rgb "#228B22";\
	set title "Quadratic spline derivatives";\
	set ylabel "y prime";\
	plot "$<" using 1:5 with lines notitle linecolor rgb "#008080",\
	"$<" using 1:6 with lines notitle linecolor rgb "#DC143C",\
	"$<" using 1:7 with lines notitle linecolor rgb "#228B22";\
	set title "Quadratic spline integrals";\
	set ylabel "Y";\
	plot "$<" using 1:8 with lines notitle linecolor rgb "#008080",\
	"$<" using 1:9 with lines notitle linecolor rgb "#DC143C",\
	"$<" using 1:10 with lines notitle linecolor rgb "#228B22";\
	unset multiplot;\
	' | tee log.gpi | gnuplot

Out.Cubic_interp.svg: out_C_data_cubic_1.txt out_C_data_cubic_2.txt Makefile
	echo '\
                set terminal svg background "white" ;\
                set out "$@" ;\
                set xlabel "x";\
                set ylabel "y";\
		set size 1,1;\
		set key outside;\
		set term svg size 600, 800 ;\
                set multiplot layout 3,1 columns ;\
                set title "Cubic spline";\
                plot "out_C_data_cubic_2.txt" with p pt 7 ps 0.5 title "Data"\
		,"" smooth csplines t "gnuplot spline"\
                ,"$<" using 1:2 with l title "Homemade spline";\
                set title "Cubic spline derivative";\
                plot "$<" using 1:3 with l notitle;\
		set title "Cubic spline integral";\
		plot "$<" using 1:4 with l notitle;\
                unset multiplot;\
                '|tee log.gpi | gnuplot

out_A_data_linear_1.txt: main_A.exe
	mono $< 1> $@ 2> out_A_data_linear_2.txt

out_B_data_quad.txt: main_B.exe
	mono $< -output:dataquad.txt 2> $@ 1> OutqTest.txt


out_C_data_cubic_1.txt: main_C.exe
	mono $< 1> $@ 2> out_C_data_cubic_2.txt

main_A.exe: main_A.cs spline.dll extlib.dll; $(MKEXE)

main_B.exe: main_B.cs spline.dll extlib.dll; $(MKEXE)

main_C.exe: main_C.cs spline.dll extlib.dll; $(MKEXE)

spline.dll: spline.cs extlib.dll; $(MKLIB)
	
extlib.dll: ../lib/matrix/matrix.cs ../lib/matrix/vector.cs\
       	../lib/LinEq/QRGS.cs ../lib/LinEq/linsol.cs \
	../lib/genlist/genlist.cs; $(MKLIB)

clean: 
	$(RM) *.dll *.exe [Oo]ut* [Ll]og*
