using System;
using static System.Console;
using static System.Math; 

class epsilon {
	public static void Main(){

		// Finding max int
		int i=1;
		while(i+1>1){i++;}
		Write($"The max int = {i}\n");
		// Finding min int
		i=1;
		while(i-1<i){i--;}
		Write($"The min int = {i}\n\n");

		//Machine Epsilon
		Write($"Machine epsilon:\n");
		double x=1;
		while(1+x!=1){x/=2;}
		x*=2;
		Write($"double= {x}\t");
		Write($"should equal Pow(2,-52)={Pow(2,-52)}\n");
		
		float y=1F; 
		while((float)(1F+y)!=1F){y/=2F;}
		y*=2F;
		Write($"float={y}\t");
		Write($"should equal Pow(2,-23)={Pow(2,-23)}\n\n");
		
		//Tinyepsilon
		Write($"Tiny Epsilon:\n");
		
		double epsilon = Pow(2,-52);
		double tiny=epsilon/2;
		double a=1+tiny+tiny;
		double b=tiny+tiny+1;
		Write($"a==b ? {a==b}\n");
                Write($"a>1 ? {a>1}\n");
                Write($"b>1 ? {b>1}\n");

		Write($"a==1 {a==1}\t b==1 {b==1}\n");
		Write("The reason a!=b is that the order of operations matter due to the rounding. Adding tiny after one causes it to round down, while adding the 2 tiny first causes the pc to round up and causes b to become slightly bigger than 1\n\n");

		Write("Comparing doubles:\n");
		double d1 = 0.1+0.1+0.1+0.1+0.1+0.1+0.1+0.1;		
		double d2 = 0.1*8;
		WriteLine($"d1={d1:e15}");
		WriteLine($"d2={d2:e15}");
		WriteLine($"d1==d2 ? => {d1==d2}");


		bool compare = approx(d1,d2);
		Write($"d1 ~ d2 : {compare}\n");
		}


	


	public static bool approx(double a, double b, double acc=1e-9, double eps=1e-9){		
		if(Abs(a-b)<=acc){return true;}
		if(Abs(a-b)<=eps*Max(Abs(a),Abs(b))){return true;}
		return false;
	}

}
