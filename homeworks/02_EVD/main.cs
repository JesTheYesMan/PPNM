using static System.Console;

public class main {
    public static void Main() {
        int matrix_size = 3;
        // Create and initialize A
        matrix matrix_A = new matrix(matrix_size);
        System.Random rng = new System.Random();
        
        for (int row = 0; row < matrix_size; row++) {
            for (int col = row; col < matrix_size; col++) {
                double random_value = rng.NextDouble();
                matrix_A[row, col] = random_value;
                matrix_A[col, row] = random_value;
            }
        }
        
        // Print matrix A
        matrix_A.print("A = ");
        
        // Compute eigenvalues and eigenvectors
        (vector eigenvalues, matrix eigenvectors) = jacobi.cyclic(matrix_A);
        (vector eigenvalues_opt, matrix eigenvectors_opt) = jacobi.cyclic_opt(matrix_A);
        
        // Print eigenvalues and eigenvectors
        eigenvalues.print("EigenValues: ");
        eigenvalues_opt.print("Eigenvalues using opt: ");
        eigenvectors.print("Eigenvectors = ");
        eigenvectors_opt.print("Eigenvectors using opt = ");
        
        // Compute matrices for validation
        matrix diagonal_matrix = matrix.diag(eigenvalues);
        matrix transposed_eigenvectors = eigenvectors.transpose();
        matrix vtav = transposed_eigenvectors * matrix_A * eigenvectors;
        matrix vdvt = eigenvectors * diagonal_matrix * transposed_eigenvectors;
        matrix vtv = transposed_eigenvectors * eigenvectors;
        
        // Print results of matrix computations
        vtav.print("V^T*A*V = ");
        vdvt.print("V*D*VT = ");
        vtv.print("V^T*V = ");
    }
}
