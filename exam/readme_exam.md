### Overview of completed homeworks

Exam evaluation: 
6-7 points


Part A: 


This part is based upon our previous work with EVD and uses the Jacobi eigenvalue algorithm for matrix diagonalization developed during handing nr 2. as well as the matrix dll

I made a couple of helper functions for inversing and making the square root of a matrix, for the sake of keeping it simple, as well as a helper function to have S be a matrix. While this could have been done without doing so, i felt it made it easier to keep an overview. 
I also introduced a matrix C to generate matrix B to make sure it was positive definete


Part B:

Fist thing i did was to combine the formulas given and check that they were correct: 
Hij = ∫0∞dr φi(r) (-1/2 d²/dr² - 1/r) φj(r) = -3/2 Sqrt(π ​αi​βj​(αi​+βj​))^(−5/2) + 1/2​(αi​+βj​)^−1,
Nij = ∫0∞dr φi(r) φj(r) = ∫0∞dr r Exp(-αr²) r  Exp(-βr²)= 1/4 Sqrt(π) (αi+βj)^-3/2.

(α+β)

I run into some issues with reachin NaN values when minimizing and amm working on implementing failsafes 

