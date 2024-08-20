using static System.Math;

public static class QRGS {

    public static (matrix Q, matrix R) decomp(matrix A) { 
        // Returns the QR decomposition in the format (Q, R)
        int matrix_size = A.size2;
        matrix matrix_Q = A.copy();
        matrix matrix_R = new matrix(matrix_size, matrix_size);
        
        for (int i = 0; i < matrix_size; i++) {
            // Compute norm for the i-th column
            matrix_R[i, i] = matrix.norm(matrix_Q[i]);
            // Normalize the i-th column of Q
            matrix_Q[i] /= matrix_R[i, i];
            
            // Orthogonalize the subsequent columns
            for (int j = i + 1; j < matrix_size; j++) {
                matrix_R[i, j] = matrix_Q[i].dot(matrix_Q[j]);
                matrix_Q[j] -= matrix_Q[i] * matrix_R[i, j];
            }
        }
        return (matrix_Q, matrix_R);
    }

    public static vector backsub(matrix A, vector b) { 
        // Perform back substitution
        for (int i = b.size - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < b.size; j++) {
                sum += A[i, j] * b[j];
            }
            b[i] = (b[i] - sum) / A[i, i];
        }
        return b;
    }

    public static vector solve(matrix A, vector b) { 
        // Solve the linear system Ax = b using QR decomposition
        (matrix Q, matrix R) = decomp(A);
        vector y = Q.transpose() * b;
        return backsub(R, y);
    }

    public static double det(matrix A) { 
        // Compute the determinant of matrix A
        if (A.size1 == A.size2) {
            matrix matrix_R = decomp(A).Item2;
            double determinant = 1;
            for (int i = 0; i < matrix_R.size1; i++) {
                determinant *= matrix_R[i, i];
            }
            return determinant;
        }
        throw new System.ArgumentException($"det: Can't compute the determinant of a non-square matrix with size ({A.size1}, {A.size2}).");
    }

    public static matrix inv(matrix A) {
        // Compute the inverse of matrix A
        if (A.size1 == A.size2) {
            int matrix_size = A.size1;
            matrix matrix_A_inv = new matrix(matrix_size);
            (matrix Q_transposed, matrix R) = decomp(A);
            matrix Q_transposed_transposed = Q_transposed.transpose();
            
            for (int i = 0; i < matrix_size; i++) {
                matrix_A_inv[i] = backsub(R, Q_transposed_transposed[i]);
            }
            return matrix_A_inv;
        }
        throw new System.ArgumentException($"inv: Can't invert a non-square matrix with size: ({A.size1}, {A.size2})");
    }
}
