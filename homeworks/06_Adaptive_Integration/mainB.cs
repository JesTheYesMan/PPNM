using System;
using static System.Console;
using static System.Math;

public class main{
	public static int Main(){
		Func<double,double> f1 = x => Pow(x, -0.5);
		Func<double,double> f2 = x => Log(x)/Sqrt(x);
		(double Res1,double err1,int eval1) = integrate.transint(f1,0,1,1e-5);
                (double Res2,double err2,int eval2) = integrate.transint(f2,0,1,1e-5);
                
                WriteLine("Functions numerically integrated:");
		WriteLine($"			1/sqrt(x):			ln(x)/sqrt(x):  \n");
		WriteLine($"Values found:		{Res1}		{Res2}	\n");
		WriteLine("Expected values:	2				-4\n");
		WriteLine($"Errors:     		{err1}		{err2}		\n");
		WriteLine($"Old evaluations:	8580				8688	");
		WriteLine($"New Evaluations:	{eval1}				{eval2}	");
                WriteLine($"Python evaluations:	231				315");
                
                
	        return 0;	
	}//Main
}//main
