using System;
using static System.Console;
using static System.Math;

public class main{
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
		int N = 5, M = 100;
                vector x = new vector(N), y1 = new vector(N), y2 = new vector(N), y3 = new vector(N), z = new double[M];
                for(int i = 0; i < N; i++){
                        x[i] = i+1;
                        y1[i] = 1;
			y2[i] = x[i];
			y3[i] = x[i]*x[i];
                        outstream.WriteLine($"{x[i]} {y1[i]} {y2[i]} {y3[i]}");} 
		interp.qspline qs1 = new interp.qspline(x,y1); //compute spline
                interp.qspline qs2 = new interp.qspline(x,y2);
                interp.qspline qs3 = new interp.qspline(x,y3);
                for(int i = 0; i < M; i++){
                        z[i] = x[0] + i*(x[N-1]-x[0])/(M-1);
                        Error.WriteLine($"{z[i]} {qs1.evaluate(z[i])} {qs2.evaluate(z[i])} {qs3.evaluate(z[i])} {qs1.derivative(z[i])} {qs2.derivative(z[i])} {qs3.derivative(z[i])} {qs1.integral(z[i])} {qs2.integral(z[i])} {qs3.integral(z[i])}");}
                vector c1 = new vector("0 0 0 0");
		vector c2 = new vector("0 0 0 0");
		vector c3 = new vector("1 1 1 1");//analytic c's
		vector b1 = new vector("0 0 0 0");
                vector b2 = new vector("1 1 1 1");
                vector b3 = new vector("2 4 6 8");//analytic b's
		if(c1.approx(qs1.c) && c2.approx(qs2.c) && c3.approx(qs3.c) 
				&& b1.approx(qs1.b) && b2.approx(qs2.b) && b3.approx(qs3.b)){
			WriteLine("Comparison test with analytical values: Success");} 
			else {WriteLine("Comparison test with analytical values: Failure");}
		qs1.c.print("c1 spline:");
		qs1.b.print("b1 spline:");
		c1.print("c1 manual:");
                b1.print("b1 manual:");
		qs2.c.print("c2 spline:");
                qs2.b.print("b2 spline:");
                c2.print("c2 manual:");
                b2.print("b2 manual:");
		qs3.c.print("c3 spline:");
                qs3.b.print("b3 spline:");
                c3.print("c3 manual:");
                b3.print("b3 manual:");
        	outstream.Close();
		return 0;
	}//Main
}//main
