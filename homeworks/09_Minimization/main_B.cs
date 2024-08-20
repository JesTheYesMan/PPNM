using System;
using static System.Console;
using static System.Math;

public class main{
	public static int Main(){
		//read the data file
		var energy = new genlist<double>();
		var signal = new genlist<double>();
		var error  = new genlist<double>();
		var separators = new char[] {' ','\t'};
		var options = StringSplitOptions.RemoveEmptyEntries;
		do{
        		string line=Console.In.ReadLine();
        		if(line==null)break;
        		string[] words=line.Split(separators,options);
        		energy.add(double.Parse(words[0]));
        		signal.add(double.Parse(words[1]));
        		error.add(double.Parse(words[2]));
		}while(true);
		int n = energy.size;
		
		Func<double,double,double,double,double> F = delegate(double E, double m, double gam, double A){
						return A/(Pow(E-m,2) + gam*gam/4);};

		vector geuss = new vector("127 3 11");
		Func<vector,double> chi2 = delegate(vector v){double m = v[0], gam = v[1], A = v[2], chi = 0; 
						for(int i=0;i<n;i++)chi+=Pow((F(energy[i],m,gam,A)-signal[i])/error[i],2);
						return chi;};
		
		min.newton minimum = new min.newton(chi2, geuss);
		vector param = minimum.x;
		double chival = minimum.f;
		int steps = minimum.steps, f_eval = minimum.f_eval;
		int N = 1000;
		double En = 0, fit = 0;
		for(int i=0;i<N;i++){
			En = energy[0] + (energy[n-1]-energy[0])*i/N;
			fit = F(En,param[0],param[1],param[2]);
			Out.WriteLine($"{En} {fit}");
		}
		Error.WriteLine($"Chi^2 value: {chival}. Number of steps: {steps}. Number of evaluations: {f_eval}");
		Error.WriteLine($"Parameters: m = {param[0]}, gamma = {param[1]}, A = {param[2]}");
		return 0;
	}//Main
}//main
