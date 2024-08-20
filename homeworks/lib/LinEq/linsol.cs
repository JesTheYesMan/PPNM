using static System.Math;

public static class linsol{

	public static vector TriDiagSol(matrix A, vector b){ //solves a tri-diagonal system in-place
		if(A.size1 != A.size2 || A.size2 != b.size){
			throw new System.ArgumentException($"TriDiagSol: Bad dimensions, A: ({A.size1},{A.size2}), b: {b.size}");}
		int n = b.size;
		for(int i = 1; i < n; i++){ //Do gaussian elimination
			double w = A[i,i-1]/A[i-1,i-1]; //Ai = A[i,i-1], Di = A[i,i], Qi = A[i,i+1]
			A[i,i] -= w*A[i-1,i];
			b[i] -= w*b[i-1];}
		vector sol = new vector(n); //Construct solution
		sol[n-1] = b[n-1]/A[n-1,n-1];
		for(int i = n-2; i >= 0; i--)sol[i] = (b[i] - A[i,i+1]*sol[i+1])/A[i,i];
		return sol;
	}//TriDiagSol
}//linsol
