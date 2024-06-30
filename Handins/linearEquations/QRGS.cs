using static System.Math;
using static System.Console;
public static class QRGS{	
	

	public static vector[] VecFromMatrix(matrix A){
		vector[] columns = new vector[A.size2];
		for (int i = 0; i < A.size2; i++){
			vector e = new vector(A.size1);
			for (int j = 0; j < A.size1; j++){
				e[j] = A[j, i];
			}
			columns[i] = e;
		}
		return columns;
	}
	
	public static (matrix, matrix) decomp(matrix A){ 
		int n = A.size1;
		int m = A.size2;
		matrix Q=new matrix(n,m), R=new matrix(m,m);
		Q=A;
		for(int i = 0; i<m; i++){
			vector ai =VecFromMatrix(A)[i]; 
			R[i,i]=ai.norm();
			Q[i]=ai/R[i,i]; //normalize the Q vectors
			for(int j=i+1; j<m; j++){
				R[i,j]=Q[i].dot(Q[j]);
				Q[j] =Q[j]-( Q[i]*R[i,j]);
			}
		}
		return (Q, R);
	}	


	public static vector solve(matrix A, vector b){ 
		(matrix Q, matrix R) = decomp(A);
		vector QTb = Q.transpose()*b;

		for(int i =  QTb.size-1; i >= 0; i--){
                        double sum = 0;
                        for(int j = i + 1; j <  QTb.size; j++)sum+=R[i,j]* QTb[j];
                        QTb[i] = (QTb[i] - sum)/R[i,i];
			}
                return QTb;
	}
	


	public static double det(matrix A){
		if(A.size1!=A.size2){
			throw new System.ArgumentException($"det: Can't take determinant of non-square matrix with size ({A.size1}, {A.size2}).");
		}
		matrix R = decomp(A).Item2;
		int size = R.size1;
		double det = 1;
		for (int i = 0; i<size; i++){
			det *= R[i, i];
		}
		return det;
	}


	public static matrix VecToMatrix(vector[] A){
	
		int n = A[0].size;
		int m = A.Length;
		matrix res = new matrix(n,m);
		for (int i=0; i<m; i++){
			for (int j=0; j<n; j++){
				res[j, i] = A[i][j];
			}
		}
		return res;
	}
	public static matrix inverse(matrix A){
		if(A.size1!=A.size2){
			throw new System.ArgumentException($"det: Can't inverse non-square matrix with size ({A.size1}, {A.size2}).");
		}

		(matrix Q, matrix R) = decomp(A);
		matrix B = new matrix(R.size1, R.size2);
		matrix QT = Q.T;
		vector[] xs = new vector[R.size1];

		for (int i=0; i<R.size1; i++){
			vector E = new vector(R.size2);
			E[i] = 1;
			xs[i] = solve(A, E);
			}

		B = VecToMatrix(xs);

		return B;
	}
}





