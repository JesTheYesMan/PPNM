using static System.Math;

public static class QRGS{
	
	public static (matrix, matrix) decomp(matrix A){ //returns in format [Q, R]
		int m = A.size2;
		matrix Q=A.copy(), R=new matrix(m,m);
		for(int i = 0; i<m; i++){
			R[i,i]=matrix.norm(Q[i]);
			Q[i]/=R[i,i]; //normalize the Q vectors
			for(int j=i+1; j<m; j++){
				R[i,j]=Q[i].dot(Q[j]);
				Q[j] -= Q[i]*R[i,j]; //Orthognalize the Q vectors
			}}
		return (Q, R);
	}//decomp	
	
	public static vector backsub(matrix A, vector b){ //back substitution inplace
		for(int i = b.size-1; i >= 0; i--){
                        double sum = 0;
                        for(int j = i + 1; j < b.size; j++)sum+=A[i,j]*b[j];
                        b[i] = (b[i] - sum)/A[i,i];
                	}
                return b;
	}//backsub

	public static vector solve(matrix A, vector b){ 
		(matrix Q, matrix R) = decomp(A);
		vector sol = Q.transpose()*b;
		return backsub(R, sol);
	}//solve
	
	public static double det(matrix A){ //only returns determinant up to sign
		if(A.size1 == A.size2){
			matrix R = decomp(A).Item2;
			double res = 1;
			for(int i = 0; i < R.size1; i++)res*=R[i,i];
			return res;
		}
		throw new System.ArgumentException($"det: Can't take determinant of non-square matrix with size ({A.size1}, {A.size2}).");
	}//det
	 
	public static matrix inv(matrix A){
		if(A.size1 == A.size2){
			int n = A.size1;
			matrix Ainv = new matrix(n);
			(matrix QT, matrix R) = decomp(A);
			QT = QT.transpose();
			for(int i = 0; i < n; i++)Ainv[i] = backsub(R,QT[i]);
		       	return Ainv;	
		}
		throw new System.ArgumentException($"inv: Can't invert non square matrix with size: ({A.size1}, {A.size2})");
	}//inv
}//QRGS
