using System;
using static System.Console;
using static System.Math;

public static class main{
	public static void Main(){
		Func<double,vector,vector> Fharm = delegate(double x, vector v){vector w = new vector(2); 
			w[0] = v[1]; w[1] = -v[0]; return w;};
		double b = 0.25, c = 5.0;
		Func<double,vector,vector> Fpend = delegate(double t, vector theta){vector w = new vector(2);
			w[0] = theta[1]; w[1] = -b*theta[1] - c*Sin(theta[0]); return w;}; //theta = (theta, omega)
		double x0 = 0, end = 10; vector u0 = new vector("3.04 0");
		var xs = new genlist<double>(); 
		var ts = new genlist<double>();
		var us = new genlist<vector>();
	        var thetas = new genlist<vector>();
		ODE.driver(Fharm, x0, u0, end, xlist: xs, ylist: us);
		ODE.driver(Fpend, x0, u0, end, xlist: ts, ylist: thetas);
		for(int i = 0; i < xs.size; i++)Out.WriteLine($"{xs[i]} {us[i][0]} ");
		for(int i = 0; i < ts.size; i++)Error.WriteLine($"{ts[i]} {thetas[i][0]} {thetas[i][1]}");
	}//Main
}//main
