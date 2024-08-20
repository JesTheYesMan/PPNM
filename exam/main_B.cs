using static System.Console;
using static System.Math;
using System;


public class main_B{
    
    	static (matrix N, matrix H) SetupMatrices(vector alpha)
    	{
        int n = alpha.size;
        int matrix_size=n;
        matrix H = new matrix(matrix_size);
        matrix N = new matrix(matrix_size);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                N[i, j] = 1.0/4.0 *Sqrt(PI)*Pow((alpha[i]+alpha[j]),-3.0/2.0);
                if (N[i, j] < 1e-10) { N[i, i] = 1e-10;}
                H[i, j] = -3.0/2.0 *Pow(Sqrt(PI*alpha[i]*alpha[j]*(alpha[i]+alpha[j])),-5.0/2.0) + 1.0/2.0*Pow((alpha[i]+alpha[j]),-1.0);
                if (H[i, j] < 1e-10) { N[i, i] = 1e-10;}
            }
        }
//    	N.print();
//    	H.print();
    	return (N,H);
    	}
    
    
    
	public static double ObjectiveFunction(vector alpha)
	{
    	(matrix N, matrix H) = SetupMatrices(alpha);
    	var (eigenvalues, _) = main_A.generalised_eigenvalues_solver(H, N);
    	//double penalty = 0.01 * alpha.dot(alpha); // Regularization term
    	//WriteLine(eigenvalues[0]);
    	return eigenvalues[0] ;//+ penalty;
	}
    
    

public static vector MinimizeAlpha(vector initialGuess)
{
    
    min.newton optimizer = new min.newton(ObjectiveFunction, initialGuess,  acc: 0.001, max_steps: 9999);
    
    
    for (int i = 0; i < optimizer.x.size; i++)
    {
        optimizer.x[i] = Max(0.1, Min(optimizer.x[i], 10.0)); //stabiliser
    }
    
    vector optimalAlphas = optimizer.x;
    double optimalEnergy = optimizer.f;
    Console.WriteLine("Optimal Alphas: ");
    optimalAlphas.print();
    Console.WriteLine("Optimisation Ground State Energy: " + optimalEnergy);
    
    return optimalAlphas;
}



    public static void Main()
    {
        
        vector initialGuess = new vector(1.5, 2.5, 3.5); // Initial guess
        initialGuess.print("Initial alpha guess");
        vector alpha = MinimizeAlpha(initialGuess );  // Example choice of variational parameters
        (matrix H, matrix N) = SetupMatrices(alpha);
        
        (vector epsilon, matrix c)=main_A.generalised_eigenvalues_solver(H, N);
        
        WriteLine("Ground State Energy: " + epsilon[0]);
        WriteLine("Ground State Vector: " );
        c[0].print();
    }
}


