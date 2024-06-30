SDKDIR=$(HOME)/.dotnet/sdk/6.0.102# or whatever your sdk is
CSCPATH=$(SDKDIR)/Roslyn/bincore/csc.dll
NETSTANDARD=$(SDKDIR)/ref/netstandard.dll
CONFIG=$(SDKDIR)/vstest.console.runtimeconfig.json
DOTNET=DOTNET_CLI_TELEMETRY_OPTOUT=1; dotnet
CSC=$(DOTNET) $(CSCPATH) -reference:$(NETSTANDARD)
RUN=$(DOTNET) exec --runtimeconfig $(CONFIG)
CSC = mcs
RUN = mono
CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) # macros
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)
TIME = time --output=$@ --append --format "$$N %e %U"



all: Out.txt; cat $<

default: Out.txt Outtimes.svg Makefile

Out.txt: main.exe
	$(RUN) main.exe -reference:matrix.dll -reference:QRGS.dll > $@

main.exe: main.cs matrix.dll QRGS.dll
	$(CSC) -reference:matrix.dll -reference:QRGS.dll  $<

matrix.dll: vector.cs matrix.cs
	$(CSC) -target:library -out:$@ $^

vector.exe: vector.cs
	$(RUN) $< > $@


OutTimes.svg: outtimes.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "Matrix size (N)";\
		set ylabel "Operational time O(N)";\
		set title "Running time depending on matreix size";\
		f(x)=b+(x/a)**c;\
		b=1; a=100; c=3;\
		fit f(x) "$<" via a,b;\
		plot "$<" with lp \
		,f(x) title sprintf("fit: (n/%.f)^{%.2f}+%.3f",a,c,b) \
		'|tee log.gpi |gnuplot


out.times.data : main.exe
	>$@
	for N in $$(seq 100 20 200); do \
		time --format "$$N %e" --output $@ --append \
		mono $< -size:$$N 1>out 2>err ;\
	done

OutTimes.txt: time.exe Makefile
	>$@
	for N in $$(seq 200 25 600); do \
		$(TIME) mono $< -msize:$$N 1>out 2>err ;\
	done ; wait

time.exe: time.cs matrix.dll QRGS.dll; $(MKEXE)

QRGS.dll: QRGS.cs matrix.dll; $(MKLIB)


clean:
	$(RM) *.dll *.exe [Oo]ut*
