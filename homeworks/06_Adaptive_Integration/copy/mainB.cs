using System;
using static System.Console;
using static System.Math;

public class main{
	public static int Main(){
		Func<double,double> f1 = x => Pow(x, -0.5);
		Func<double,double> f2 = x => Log(x)/Sqrt(x);
		(double Res1,double err1,int eval1) = integrate.transint(f1,0,1,1e-5);
                (double Res2,double err2,int eval2) = integrate.transint(f2,0,1,1e-5);
		WriteLine("Python evaluations are computed using scipy.integrate.quad routine with same requested precision.");
                WriteLine("");
		WriteLine($"1: f(x) = 1/sqrt(x); Analytic: 2, Numerical: {Res1}, err: {err1}, #evals: {eval1}, python evals: 231");
                WriteLine($"2: f(x) = log(x)/sqrt(x); Analytic: -4, Numerical: {Res2}, err: {err2}, #evals: {eval2}, python evals: 315");
	        return 0;	
	}//Main
}//main
