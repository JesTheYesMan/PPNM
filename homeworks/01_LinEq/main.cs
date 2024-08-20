using static System.Console;

static class main {
    static int Main() {
        // Initialize the random number generator
        System.Random rng = new System.Random();
        
        // Define matrix dimensions
        int num_rows = 4;
        int num_cols = 3;
        
        // Create matrix A with the specified number of rows and columns
        matrix matrix_A = new matrix(num_rows, num_cols);
        
        // Fill matrix A with random values between 0 and 1
        for (int row = 0; row < num_rows; row++) {
            for (int col = 0; col < num_cols; col++) {
                matrix_A[row, col] = rng.NextDouble();
            }
        }
        
        WriteLine("Part A: Test QR factorization");

        // Print matrix A
        matrix_A.print("A = ");
        
        // Perform QR decomposition on matrix A
        (matrix q_matrix, matrix r_matrix) = QRGS.decomp(matrix_A);
        
        // Print Q and R matrices
        q_matrix.print("Q = ");
        r_matrix.print("R = ");
        
        // Calculate Q^T * Q
        matrix q_transpose_Q = q_matrix.transpose() * q_matrix;
        q_transpose_Q.print("Q^TQ =");
        
        // Calculate QR product
        matrix qr_product = q_matrix * r_matrix;
        qr_product.print("QR =");
        
        WriteLine("");
        
        // Verify QR decomposition and orthogonality of Q
        bool is_decomposition_correct = qr_product.approx(matrix_A) && q_transpose_Q.approx(matrix.id(num_cols));
        if (is_decomposition_correct) {
            WriteLine("Test; QR = A, Q orthogonal: Success");
        } else {
            WriteLine("Test; QR = A, Q orthogonal: Failure");
        }
        
        WriteLine("");
        WriteLine("Test linear equation solver on square matrix");
        
        // Create and fill matrix B and vector b with random values
        matrix matrix_B = new matrix(num_rows);
        vector vector_b = new vector(num_rows);
        for (int i = 0; i < num_rows; i++) {
            vector_b[i] = 2 * rng.NextDouble() - 1;
            for (int j = 0; j < num_rows; j++) {
                matrix_B[i, j] = 2 * rng.NextDouble() - 1;
            }
        }
        
        // Print matrix B and vector b
        matrix_B.print("B = ");
        vector_b.print("b = ");
        
        // Solve the linear system Bx = b
        vector solution = QRGS.solve(matrix_B, vector_b);
        solution.print("x = ");
        
        // Compute B * x
        vector computed_Bx = matrix_B * solution;
        computed_Bx.print("Bx = ");
        
        // Verify the solution
        if (computed_Bx.approx(vector_b)) {
            WriteLine("Test of linear equation solver: Success");
        } else {
            WriteLine("Test of linear equation solver: Failure");
        }
        
        WriteLine("");
        WriteLine("Part B: Matrix inverse:");
        WriteLine("");
        
        // Create and fill matrix C with random values
        matrix matrix_C = new matrix(num_rows);
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_rows; j++) {
                matrix_C[i, j] = rng.NextDouble();
            }
        }
        
        // Compute the inverse of matrix C
        matrix matrix_C_inverse = QRGS.inv(matrix_C);
        
        // Print matrix C and its inverse
        matrix_C.print("C = ");
        matrix_C_inverse.print("C^-1 = ");
        
        // Compute C * C^-1
        matrix identity_matrix = matrix_C * matrix_C_inverse;
        identity_matrix.print("C*C^-1 = ");
        
        // Verify the matrix inverse
        if (identity_matrix.approx(matrix.id(num_rows))) {
            WriteLine("Test of matrix inverse: Success");
        } else {
            WriteLine("Test of matrix inverse: Failure");
        }
        
        return 0;
    }
}
