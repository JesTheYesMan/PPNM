Out1.txt : main1.exe
	mono main1.exe  -numbers:1,2,3,4,5 > Out1.txt

Out2.txt : main2.exe
	echo 1 2 3 4 5 > input.txt
	mono main2.exe < input.txt 2> Out2.txt
Out3.txt : main3.exe
	mono main3.exe -input:input2.txt -output:Out3.txt


main1.exe: main1.cs
	mcs -target:exe -out:$@ $<

main2.exe: main2.cs
	mcs -target:exe -out:$@ $<

main3.exe: main3.cs
	mcs -target:exe -out:$@ $<
