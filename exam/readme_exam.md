### Overview of completed homeworks

Exam evaluation: 
6-7 points
I believe to have achieved all 6 points on the first task, but cant get the second to fully work

Part A: 


This part is based upon our previous work with EVD and uses the Jacobi eigenvalue algorithm for matrix diagonalization developed during handing nr 2. as well as the matrix dll

I made a couple of helper functions for inversing and making the square root of a matrix, for the sake of keeping it simple, as well as a helper function to have S be a matrix. While this could have been done without doing so, i felt it made it easier to keep an overview. 
I also introduced a matrix C to generate matrix B to make sure it was positive definete


Part B:

Fist thing i did was to combine the formulas given and check that they were correct: 
Hij = ∫0∞dr φi(r) (-1/2 d²/dr² - 1/r) φj(r) = -3/2 Sqrt(π ​αi​βj​(αi​+βj​))^(−5/2) + 1/2​(αi​+βj​)^−1,
Nij = ∫0∞dr φi(r) φj(r) = ∫0∞dr r Exp(-αr²) r  Exp(-βr²)= 1/4 Sqrt(π) (αi+βj)^-3/2.

I turn these into a function that can create these matrices from a given alpha 
I then create an objective function and a minimize function

I run into some issues with reaching NaN values in the new H and N when minimizing and tried implementing failsafes 
I have tried making alot of restriction on some of the values such that the different falues dont become too low or high an cause instability, but that hasnt worked
When trying to use qnewton the minimisation for some reason gets stuck, so can only use newton minimisation which isnt ideal. This means that some
This is unfortunately where i got stuck and i cant seam to get past this hurdle. 


with the starting guess do i get: 

Initial alpha guess       1.5        2.5        3.5 
Optimal Alphas: 
      1.83       3.09       4.35 
c: -0.210550752245064
Ground State Energy: 0.216708594386718
Ground State Vector: 
      4.81      -23.8       24.6 
      
This result seems unreliable however

I assume either something is wrong with the way the calculation of the Optimisation Ground State Energy or that the instability doesnt allow the minimization to converge. 


