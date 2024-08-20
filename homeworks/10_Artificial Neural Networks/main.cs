using System;
using static System.Console;
using static System.Math;

public static class main{

	public static (vector, vector) bound(int n, vector min, vector max){
		vector pmin = new vector(3*n);
		vector pmax = new vector(3*n);
		for(int i=0;i<n;i++){
			pmin[3*i+0] = min[0]; pmax[3*i+0] = max[0];
                        pmin[3*i+1] = min[1]; pmax[3*i+1] = max[1];
                        pmin[3*i+2] = min[2]; pmax[3*i+2] = max[2];
		}
		return (pmin,pmax);
	}

	public static int Main(string[] args){
		//Config extra output file
		string outfile=null;
		foreach(var arg in args){
			string[] words = arg.Split(':');
			if(words[0]=="-output")outfile=words[1];
		}
		if(outfile==null){
			System.Console.Error.WriteLine("wrong argument");
			return 1;
		}
		var outstream = new System.IO.StreamWriter(outfile,append:false);
		
		//start of actual program
		int N = 30; /* number of sample points */	
		int nodes = 3;
		vector max = new vector("10 10 10"), min = new vector("-10 0.001 -10");	// parameter bounds for global optimization
		Func<double,double> g = X => Sin(4*X-1)*Exp(-X*X);
		vector x = new vector(N), y = new vector(N);
		for(int i=0;i<N;i++){
			x[i] = -1+2.0*i/(N-1); y[i] = g(x[i]);
			outstream.WriteLine($"{x[i]} {y[i]}");}
		ann.interpol[] Networks = new ann.interpol[3];
		for(int i=0;i<3;i++){
			Networks[i] = new ann.interpol(nodes);
			(vector pmin, vector pmax) = bound(nodes,max,min);
			Networks[i].train_interp(x,y,pmin:pmin,pmax:pmax);
			Error.WriteLine($"Training with {nodes} nodes terminated with success: {Networks[i].train_status}");
		}
		int M = 100;
		double z = -1;
		for(int i=0;i<M;i++){
			Out.WriteLine($"{z} {Networks[0].response(z)} {Networks[0].dresponse(z)} {Networks[0].Iresponse(-1,z)} {Networks[1].response(z)} {Networks[1].dresponse(z)} {Networks[1].Iresponse(-1,z)} {Networks[2].response(z)} {Networks[2].dresponse(z)} {Networks[2].Iresponse(-1,z)}");
			z+=2.0/(M-1);
		}
		outstream.Close();
		return 0;
	}//Main
}//main
