using static System.Math;

namespace matrix_decomp{

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

public class LU{
	
	public readonly int dim;
	public readonly matrix L, U, A;
	
	//constructor
	
	public LU(matrix B){
		if(B.size1!=B.size2) throw new System.ArgumentException($"LU: non-square matrix A: ({B.size1},{B.size2})");
		A = B.copy();
		dim = A.size1;
		L = matrix.id(dim); U = new matrix(dim);
		for(int i=0;i<dim;i++){
			for(int j=i;j<dim;j++){
				double sum = 0;
				for(int k=0;k<j;k++)sum+=L[i,k]*U[k,j];
				U[i,j]=A[i,j] - sum;
			}	
			for(int j = i+1;j<dim;j++){
                                double sum = 0;
                                for(int k=0;k<j;k++)sum+=L[j,k]*U[k,i];
				L[j,i] = (A[j,i]-sum)/U[i,i];
			}	
		}
	}
	
	public double determinant(){
                double res = 1;
                for(int i=0;i<dim;i++)res*=U[i,i];
                return res;
                }

	vector backsub(vector b){ //back substitution inplace
        	for(int i = dim-1; i >= 0; i--){
                	double sum = 0;
                        for(int j = i + 1; j < dim; j++)sum+=U[i,j]*b[j];
                        	b[i] = (b[i] - sum)/U[i,i];
                        }
                return b;
                }//backsub

         vector forwardsub(vector b){ //forward substitution inplace
                for(int i = 0; i < dim; i++){
                        double sum = 0;
                        for(int k = 0; k < i; k++)sum+=L[i,k]*b[k];
         	               b[i] = (b[i] - sum)/L[i,i];
                        }
                return b;
                }//backsub

         public vector linsol(vector b){ //solves the problem inplace
         	return backsub(forwardsub(b));
                }

         public matrix inverse(){
         	matrix Ainv = matrix.id(dim);
                for(int i=0;i<dim;i++)linsol(Ainv[i]);
                return Ainv;
                }
	
	}//LU
public class cholesky{
		public readonly matrix A,L,LT;
		public readonly int dim;

		//constructor
		public cholesky(matrix B){
			if(B.size1 != B.size2) throw new System.ArgumentException($"LLT: non-square matrix, size: ({B.size1}, {B.size2})");
			dim = B.size1;
			L = new matrix(dim); A = B;
			for(int i=0;i<dim;i++)
			for(int j=0;j<=i;j++){
				double sum = 0;
				for(int k=0;k<j;k++)sum+=L[i,k]*L[j,k];
				if(i==j) L[i,j] = Pow(B[i,j]-sum,0.5);
				else L[i,j] = (B[i,j]-sum)/L[j,j];
			}
			LT = L.transpose();
		}

		//methods

		public double determinant(){
			double res = 1;
			for(int i=0;i<dim;i++)res*=L[i,i];
			return res*res;
		}

		vector backsub(vector b){ //back substitution inplace
			for(int i = dim-1; i >= 0; i--){
	                        double sum = 0;
        	                for(int j = i + 1; j < dim; j++)sum+=LT[i,j]*b[j];
	                        b[i] = (b[i] - sum)/LT[i,i];
                		}
                	return b;
		}//backsub

		vector forwardsub(vector b){ //forward substitution inplace
                        for(int i = 0; i < dim; i++){
                                double sum = 0;
                                for(int k = 0; k < i; k++)sum+=L[i,k]*b[k];
                                b[i] = (b[i] - sum)/L[i,i];
                                }
                        return b;
                }//backsub
		
		public vector linsol(vector b){ //solves the problem inplace
			return backsub(forwardsub(b));
		}

		public matrix inverse(){
			matrix Ainv = matrix.id(dim);
			for(int i=0;i<dim;i++)linsol(Ainv[i]);
			return Ainv;
		}

	}//cholesky
}//matrix decomp
