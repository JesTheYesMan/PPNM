CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)
TIME = time --output=$@ --append --format "$$N %e %U"

all: Out.txt Out.fit.svg

Out.fit.svg: data.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@";\
		set xlabel "t (days)" ;\
		set ylabel "Decays" ;\
		set grid;\
		set title "exponential fit of radioactive decay of 224Ra";\
		f(x) = a*exp(-b*x) ;\
		fp(x) = a*exp(da)*exp(-(b-db)*x) ;\
		fm(x) = a*exp(-da)*exp(-(b+db)*x) ;\
		a = 142.403540630783; b = 0.171; da = 0.0354019292641024; db = 0.0071566223420022;\
		plot "$<" with errorbars title "data"\
		,f(x) title sprintf("fit: A(t)= %.f*e^{-%.2f t}",a,b)\
		,fp(x) with lines title sprintf("fit with uncertainty upper bound ")\
		,fm(x) with lines title sprintf("fit with uncertainty lower bound")\
		'|tee log.gpi | gnuplot

Out.txt: main.exe
	mono $< -input:data.txt > $@

main.exe: main.cs extlib.dll lsfit.dll; $(MKEXE)

extlib.dll: ../lib/matrix/matrix.cs ../lib/matrix/vector.cs \
	../lib/LinEq/QRGS.cs ../lib/genlist/genlist.cs; $(MKLIB) 
	
lsfit.dll: extlib.dll lsfit.cs; $(MKLIB)

clean: 
	$(RM) *.dll *.exe [Oo]ut* [Ll]og*
