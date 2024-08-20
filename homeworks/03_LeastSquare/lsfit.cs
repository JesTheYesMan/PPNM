using static System.Math;

public static class fit {
    
    public static (vector popt, matrix pcov) lsfit(
        System.Func<double, double>[] functions, 
        vector x_values, 
        vector y_values, 
        vector uncertainties) {
        
        // Check if input data sizes are consistent
        if (x_values.size != y_values.size || x_values.size != uncertainties.size || y_values.size != uncertainties.size) {
            throw new System.ArgumentException($"lsfit: Incompatible data sizes: ({x_values.size} {y_values.size} {uncertainties.size}).");
        }
        
        // Initialize matrices and vectors
        matrix design_matrix = new matrix(x_values.size, functions.Length);
        vector b_vector = new vector(y_values.size);
        
        // Fill the design matrix and b vector
        for (int i = 0; i < x_values.size; i++) {
            b_vector[i] = y_values[i] / uncertainties[i];
            for (int k = 0; k < functions.Length; k++) {
                design_matrix[i, k] = functions[k](x_values[i]) / uncertainties[i];
            }
        }
        
        // Perform QR decomposition
        (matrix Q, matrix R) = QRGS.decomp(design_matrix);
        
        // Solve for the optimal parameters
        vector popt = QRGS.backsub(R, Q.transpose() * b_vector);
        
        // Compute the parameter covariance matrix
        matrix R_inv = QRGS.inv(R);
        matrix pcov = R_inv * R_inv.transpose();
        
        return (popt, pcov);
    } // lsfit
} // fit
