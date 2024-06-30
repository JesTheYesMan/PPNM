using static System.Console;
using System;

class main{
	static int Main(){
		System.Random random = new System.Random();
		int m = 5, n = 3;
		matrix A = new matrix(m,n);
		for(int i = 0; i < m; i++){
			for(int j = 0; j < n; j++){ 
				A[i,j] = random.NextDouble();
				}
			}

		WriteLine("Task A Testing QR factorization");

		A.print("A = ");
		(matrix Q, matrix R) = QRGS.decomp(A);
		Q.print("Q = ");
		R.print("R = ");
		matrix QQ = Q.transpose()*Q;
		QQ.print("Q^TQ =");
		matrix QR = Q*R;
		QR.print("QR =");

		WriteLine("");
		if(QR.approx(A) && QQ.approx(matrix.id(n))) WriteLine("Test; QR = A, Q orthogonal: Success");
                else WriteLine("Test; QR = A, Q orthogonal: Failure");
		WriteLine("");
		WriteLine("Test linear equation solver on square matrix");
		matrix B = new matrix(m);
		vector b = new vector(m);
                for(int i = 0; i < m; i++){ 
			b[i] = 2*random.NextDouble()-1;
			for(int j = 0; j < m; j++) B[i,j] = 2*random.NextDouble()-1;}
		B.print("B = ");
		b.print("b = ");
		vector sol = QRGS.solve(B, b);
                sol.print("x = ");
		vector Bx = B*sol;
                Bx.print("Bx = ");

		if(Bx.approx(b)) WriteLine("Test of linear equation solver: Success");
		else WriteLine("Test of linear equation solver: Failure");

		WriteLine("");
		WriteLine("Task B, Matrix inverse:");
		WriteLine("");
		matrix C = new matrix(m);
                for(int i = 0; i < m; i++){
                for(int j = 0; j < m; j++) C[i,j] = random.NextDouble();}
		matrix Cinv = QRGS.inv(C);
		C.print("C = ");
		Cinv.print("C^-1 = ");
		matrix CCinv = C*Cinv;
		CCinv.print("C*C^-1 = ");
		if(CCinv.approx(matrix.id(m))) WriteLine("Test of matrix inverse: Success");
                else WriteLine("Test of matrix inverse: Failure");
		return 0;
	}//Main
}//main
