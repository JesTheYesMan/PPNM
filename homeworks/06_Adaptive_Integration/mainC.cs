using System;
using static System.Console;
using static System.Math;

public class main{
	public static int Main(){
		Func<double,double>[] fs = new Func<double,double>[4]{x => 1/(1+Pow(x,2)), x=> Exp(-x), x=> Exp(-x)*Pow(x,2), x=> Exp(-Pow(x,2))};
		
		string[] names = new string[4]{"1/(1+x^2)", "exp(-x)", "x^2e^-x", "e^(-x^2)"};
	      	
	      	double[] vals = new double[4]{PI, 1, 2, Sqrt(PI)};
		
		double[] a = new double[4]{
			double.NegativeInfinity, 
			0, 
			0,
			double.NegativeInfinity};
		
		double[] b = new double[4]{
			double.PositiveInfinity, 
			double.PositiveInfinity, 
			double.PositiveInfinity, 
			double.PositiveInfinity};
		
		double[] res= new double[4]{0,0,0,0};
		double[] err= new double[4]{0,0,0,0};
		double[] eval= new double[4]{0,0,0,0};
		
		int[] PyN= {90,135,165,270};

		for(int i = 0; i < 4; i++){
		(res[i],err[i],eval[i]) = integrate.integral(fs[i],a[i],b[i]);
		}	

		WriteLine("Functions numerically integrated:");
		WriteLine($"			{names[0]}		{names[1]}:		{names[2]}:		{names[3]}:  \n");
		WriteLine($"limits:			-inf, inf		0,inf			0,inf			-inf,inf \n");
		WriteLine($"Values found:		{res[0]}	{res[1]}	{res[2]}	{res[3]}	\n");
		WriteLine($"Expected values:	{vals[0]}	{vals[1]}			{vals[2]}			{vals[3]}\n");
		WriteLine($"Errors:     		{err[0]}	{err[1]}	{err[2]}	{err[3]}\n");
		WriteLine($"Evaluations:		{eval[0]}			{eval[1]}			{eval[2]}			{eval[3]}\n");
                WriteLine($"Python evaluations:	{PyN[0]}			{PyN[1]}			{PyN[2]}			{PyN[3]}");
	        return 0;	
	}//Main
}//main


