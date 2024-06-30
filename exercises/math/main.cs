using System;
using static System.Console;
using static System.Math;

class main{
	public static void Main(){
		double sqrt2=Sqrt(2.0);
		Write($"sqrt(2=) {sqrt2}\n");
		Write($"sqrt2^2 ={sqrt2*sqrt2} (should eaqual 2) \n");

		double power2to15 = Pow(2.0,1.0/5.0);
		Write($"2^1/5 is {power2to15:F6} which should equal  1.148..)\n");
		Write($"(2^1/5)^5= {Pow(power2to15,5.0)} (should eaqual 2)\n");

		double PowerEToPi = Exp(PI);
		Write($"e^pi = {PowerEToPi:F6} (should equal 23.14...)\n");

		double PowerPiToE = Pow(PI,E);
		Write($"pi^e = {PowerPiToE:F6} (should eaqual 22.45...)\n");
		
		for(int i=1;i<=10;i++){
			Write($"For {i-1}:\t\t");
			Write($"Gamma({i})= {sfuns.fgamma(i):F2}\t\t");
			Write($"ln(Gamma({i}))={sfuns.lngamma(i):F8} -> exp(ln(Gamma({i}))) = {Exp(sfuns.lngamma(i)):F8}\n");
		}
	}
}
