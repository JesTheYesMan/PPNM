using System;
using static System.Console;
using static System.Math;

static class main{
	static void Main(){
		matrix M = new matrix("1 0 0 0; 0 1 0 0; 0 0 1 0"); 
		double[] a = new double[]{1,2,3};
		double[] b = new double[]{-1, -2, -3, -4};
		matrix Ma = matrix.diag(b);
                matrix Ma1 = matrix.diag(a, 1);
                matrix Ma2 = matrix.diag(a, -1);
		Ma.print();
		Ma1.print();
		Ma2.print();
	}
}

