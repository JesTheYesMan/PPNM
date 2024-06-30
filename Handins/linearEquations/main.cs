using static vector;
using static System.Console;
public static class Matrix{
	public static matrix Random(int size1, int size2){
		matrix RandomMatrix = new matrix(size1, size2);
		var rnd = new System.Random();
		for (int i=0; i<size1; i++){
			for (int j=0; j<size2; j++){
				RandomMatrix[i, j] = rnd.NextDouble()*10 - 5;
			}
		}
		return RandomMatrix;
	}

	public static vector Random(int size){
		vector vec = new vector(size);

		var rnd = new System.Random();

		for (int i = 0; i<size; i++){
			vec[i] = rnd.NextDouble()*10-5;
		}
		return vec;
	}


public static vector[] Vecs(matrix A){
		vector[] columns = new vector[A.size2];
		for (int i = 0; i < A.size2; i++){
			vector e = new vector(A.size1);
			for (int j = 0; j < A.size1; j++){
				e[j] = A[j, i];
			}
			columns[i] = e;
		}
		return columns;
	}//Vecs

	public static matrix VecsToMatrix(vector[] A){
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

}


public static class main{
static void Main(){


	Write("Part A:\n\n");

	Write("Random matrix of size 5x3\n");
	matrix A1 = Matrix.Random(5,3);
	A1.print();

	Write("\n\n Decomp test \n");
        (matrix Q1, matrix R1) = QRGS.decomp(A1);
        Q1.print();

	Write("\n\n Checking that R is upper tringular:");
        R1.print();

	Write("\n\n Checking Q^T*Q=1:");
	matrix QTQ = Q1.T*Q1;
	QTQ.print();

	Write("\n\n Checking  QR=A\n");
	matrix QR=Q1*R1;
	QR.print();


	Write("\n\n Checking solve:\n");

	Write("Random matrix of size 3x3\n");
	matrix A = Matrix.Random(3,3);
	A.print();


	matrix I = matrix.id(A.size2);

	Write("\n\n Random vector of size 3\n");
	vector H = Matrix.Random(3);
	H.print();





	Write("\n\n Decomposing A  \n");
	(matrix Q, matrix R) = QRGS.decomp(A);
	Write("Q=");
	Q.print();
	Write("R=");
	R.print();


	matrix C = Q*R;
	Write("\n\n Checking that QR=A:\n");
	C.print();

	Write("\n Decomposing\n");
	//(matrix Q, matrix R) = QRGS.decomp(A);
	Q.print();
	Write("\n\n The solution to the system QRx=b:\n");
	vector a = QRGS.solve(A, H);
	a.print();
	Write("\n\n Checking that Ax=b: \n");
	vector D = A*a;
	D.print();

	Write("\n\n Part B:\n");

	Write("I continue using the same matrix A as before so the first 2 tasks are done:\n");
		

	//Inverse
	matrix A_inverse = QRGS.inverse(A);
	Write("calculate B");
	A_inverse.print();
	Write("\n \n Testing that B is the inverse:");
	Write($"\n A*B=A*A^-1=I : {I.approx(A*A_inverse)}");

	}
}
