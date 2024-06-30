using System;
using static System.Console;
using static System.Math;
using static matrix;

public static class jacobi{

	public static void timesJ(matrix A, int p, int q, double theta){
	double c = Math.Cos(theta), s = Math.Sin(theta);
	for(int i = 0; i < A.size1; i++){
		double aip = A[i,p], aiq = A[i,q];
		A[i,p] = c*aip - s*aiq;
		A[i,q] = s*aip + c*aiq;
		}
	}

	public static void Jtimes(matrix A, int p, int q, double theta){
	double c = Math.Cos(theta),s = Math.Sin(theta);
	for(int j = 0; j < A.size1 ; j++){
		double apj = A[p,j], aqj = A[q,j];
		A[p,j] = c*apj + s * aqj;
		A[q,j] = -s*apj + c * aqj;
		}
	}
	public static (matrix, matrix, matrix, matrix) cyclic(matrix A, matrix V){
		matrix oldA = A.copy();
		matrix oldV = V.copy();
		bool changed;
		int n = A.size1; // Assuming A is square
		do{
			changed = false;
			for (int p = 0; p < n - 1; p++){
				for (int q = p + 1; q < n; q++){
					double apq = A[p, q], app = A[p, p], aqq = A[q, q];
					double theta = 0.5 * Atan2(2 * apq, aqq - app);
					double c = Cos(theta), s = Sin(theta);
					double new_app = c * c * app - 2 * s * c * apq + s * s * aqq;
					double new_aqq = s * s * app + 2 * s * c * apq + c * c * aqq;
					if (new_app != app || new_aqq != aqq){
						changed = true;
						timesJ(A, p, q, theta); // A←A*J 
						Jtimes(A, p, q, -theta); // A←JT*A 
						timesJ(V, p, q, theta); // V←V*J
					}
				}
			}
		}
		while (changed);
		return (oldA, A, oldV, V);
	}//cyclic
}//Jacobi

public static class Matrix{
	public static matrix Random(int size1, int size2){
		matrix RandomMatrix = new matrix(size1, size2);
		var rnd = new System.Random(1);
		for (int i=0; i<size1; i++){
			for (int j=i; j<size2; j++){
				RandomMatrix[i, j] = rnd.NextDouble()*10 - 5;
				RandomMatrix[j, i] = RandomMatrix[i, j]; //To make the matrix symmetric
			}
		}
		return RandomMatrix;
	}// Random matrix
}

public static class main{
	static void Main(){
		matrix A = Matrix.Random(3,3);
		matrix v = matrix.id(A.size2);
		(matrix oldA, matrix D, matrix oldV, matrix updV) = jacobi.cyclic(A, v);
		WriteLine($"Testing A=V*D*V^T : {oldA.approx(updV*D*updV.T)}\n");
		WriteLine($"Testing D=V^T*A*V : {D.approx(updV.T*oldA*updV)}\n");
		WriteLine($"Testing V^T*V=1 : {oldV.approx(updV.T*updV)}\n");
		WriteLine($"Testing V*V^T : {oldV.approx(updV*updV.T)}\n");
	}//end of Main()
}//end of main
