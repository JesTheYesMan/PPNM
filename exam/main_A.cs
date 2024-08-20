using static System.Console;
using static System.Math;
using System;
public class main_A {

public static void Main(){

	 int matrix_size = 3;
        // Create and initialize random symmetrical A and B 
        matrix A = new matrix(matrix_size);
        matrix B = new matrix(matrix_size);
        
        
	matrix C = new matrix(matrix_size); // Matrix used to make B positive definite 
	
        System.Random rng = new System.Random();
        
        for (int row = 0; row < matrix_size; row++) {
            for (int col = row; col < matrix_size; col++) {
                double random_value = rng.NextDouble();
                A[row, col] = random_value;
                A[col, row] = random_value;
                random_value = Abs(rng.NextDouble());
                C[row, col] = random_value;
            }
        }
        
        
        B = C.transpose() * C + 1e-10 * matrix.id(matrix_size);
        
        generalised_eigenvalues_solver(A, B);
	}
	
	
public static (vector eigenvalues, matrix V) generalised_eigenvalues_solver(matrix A, matrix B) {

       	int n=A.size1;
        WriteLine("Solving generalised eigenvalue problem");
        // Print matrix A
        A.print("\nA = ");
        B.print("\nB = ");
        
        
        // Finding Q and S via diagonalisation
        var (S, Q) = Diagonalize(B);
        
        Q.print("\nQ = ");
        S.print("\nS = ");
        
        
        // Finding S^(-1/2)
    matrix S_inv = InverseDiagonal(S);
    	S_inv.print("\nS_inv = ");
    	
    matrix S_inv_sqrt = matrix_sqrt(S_inv);
        S_inv_sqrt.print("\nS_inv_sqrt = ");

	//Finding A_tilde
    matrix QTAQ = Q.transpose() * A * Q; 
    matrix A_tilde= new matrix(n);
    for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
            {
                A_tilde[i, j] = S_inv_sqrt[i,i] *S_inv_sqrt[j,j] * QTAQ[i, j];
            }
        }
    //matrix A_tilde = S_inv_sqrt * QTAQ * S_inv_sqrt;
    	A_tilde.print("\nA_tilde = ");
    
    // Solving the standard eigenvalue problem to find X and the eigenvalues
    var (eigenvalues, X) = jacobi.cyclic(A_tilde);
    
    eigenvalues.print("\nEigenvalues = ");
    
    // Finding V 
    matrix V = Q * S_inv_sqrt * X;
    
    V.print("\nV = ");
    
    return(eigenvalues, V);

}
// Helper method to diagonalize a matrix using jacobi
public static (matrix S, matrix Q) Diagonalize(matrix B)
{
    (vector eigenvalues, matrix Q) = jacobi.cyclic(B);
    matrix S= matrix.diag(eigenvalues);
    return(S,Q);
}
//Helper method to make S^(1/2)
public static matrix matrix_sqrt(matrix S)
{
	matrix sqrt_S = new matrix(S.size1);
	for(int i = 0; i < S.size1; i++){
		sqrt_S[i,i]=Sqrt(S[i,i]);
	}
	
	return sqrt_S;
}

public static matrix InverseDiagonal(matrix A)
{
    if (A.size1!=A.size2){throw new Exception("Matrix is not square");}
    int n = A.size1; 
    for(int i = 0; i < n; i++)
    for(int j = i+1; j < n; j++)
    {if (i!=j && A[i,j]!=0) {throw new Exception("Matrix is not diagonal");}}
    
    matrix inverse = new matrix(n); 

    for (int i = 0; i < n; i++)
    {
        double diagElement = A[i, i];
        if (diagElement == 0){throw new Exception("Matrix is singular and cannot be inverted.");}

        inverse[i, i] = 1.0 / diagElement; // Inverse of the diagonal element
    }

    return inverse;
}




}//mian




