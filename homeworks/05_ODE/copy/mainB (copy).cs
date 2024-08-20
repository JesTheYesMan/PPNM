using System;
using static System.Console;
using static System.Math;

class main{
	static int Main(string[] args){
		double epsilon = 0, phi0 = 0, phiL = 2*PI; vector u0 = new vector("0 0");
		foreach(var arg in args){
			string[] words = arg.Split(':');
			if(words[0]=="-epsilon")epsilon=double.Parse(words[1]);
		        if(words[0]=="-phiL") phiL = double.Parse(words[1]);
			if(words[0]=="-u0") u0[0] = double.Parse(words[1]);
			if(words[0]=="-uprime0") u0[1] = double.Parse(words[1]);
		}
		Func<double,vector,vector> f = delegate(double phi, vector u){vector w = new vector(2);
			w[0] = u[1]; w[1] = 1-u[0] + epsilon*u[0]*u[0]; return w;};
		vector ul = ODE.driver(f, phi0, u0, phiL, hmax: 10);
		interp.vcspline sol = ODE.driver_interp(f, phi0, u0, phiL, hmax: 10);
		Error.WriteLine($"ul = {ul[0]}, ul\' = {ul[1]}");
		int N = 400; double ph = phi0;
		for(int i = 0; i < N; i++){WriteLine($"{ph} {sol.evaluate(ph)[0]}"); ph+=(phiL-phi0)/N;}
		return 0;
	}//Main
}//main
