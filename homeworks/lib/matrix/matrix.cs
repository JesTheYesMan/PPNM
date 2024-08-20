using System;
using static System.Console;
using static System.Math;

public partial class matrix{

public readonly int size1, size2;
public double[][] data;

//constructors
public matrix(int n, int m){
	size1 = n; size2 = m; data = new double[size2][];
	for(int i = 0; i<size2; i++) data[i] = new double[size1];
}

public matrix(int n) {
	size1 = n; size2 = n; data = new double[n][];
	for(int i = 0; i< n; i++) data[i] = new double[n];
}

public matrix(string s){
        string[] rows = s.Split(';');
        size1 = rows.Length;
	char[] delimiters = {',',' '};
        var options = StringSplitOptions.RemoveEmptyEntries;
        size2 = rows[0].Split(delimiters,options).Length;
        data = new double[size2][];
	for(int j=0;j<size2;j++) data[j]=new double[size1];
        for(int i=0;i<size1;i++){
                string[] ws = rows[i].Split(delimiters,options);
                for(int j=0; j<size2; j++){
                        this[i,j]=double.Parse(ws[j]);
        	}
	}
}

public double this[int n, int m]{
	get{return data[m][n];}
	set{data[m][n] = value;}
}

public vector this[int n]{
	get{return (vector)data[n];}
	set{data[n] = (double[])value;}
}

public void set(int i, int j, double a){data[j][i] = a;}
public void set(int i, vector a, bool row = false){
	if(!row)data[i] = (double[])a;
	if(row)for(int j = 0; j < size2; j++)data[j][i] = a[j];	
}

public matrix transpose(){
	matrix M = new matrix(size2, size1);
	for(int i = 0; i < size2; i++){
	for(int j = 0; j < size1; j++)
		M[i,j] = this[j,i];
	}
	return M;
}

//Operators
public static matrix operator-(matrix M){
	matrix A = new matrix(M.size1, M.size2);
	for(int ir=0; ir < M.size1; ir++){
        for(int ic=0; ic < M.size2; ic++)
                A[ir, ic] = -M[ir, ic];
        }
	return A;
}

public static matrix operator+(matrix A, matrix B){
	if(A.size1 == B.size1 && A.size2 == B.size2){
        	matrix M = new matrix(A.size1, A.size2);
        	for(int ir=0; ir < M.size1; ir++){
        	for(int ic=0; ic < M.size2; ic++)
                	M[ir, ic] = A[ir, ic] + B[ir, ic];
        	}
        	return M;
		}
	throw new ArgumentException($"Matrix +: Matrices with sizes ({A.size1}, {A.size2}) and ({B.size1}, {B.size2}) cannot be added.");
}

public static matrix operator-(matrix A, matrix B){
        if(A.size1 == B.size1 && A.size2 == B.size2){
                matrix M = new matrix(A.size1, A.size2);
                for(int ir=0; ir < M.size1; ir++){
                for(int ic=0; ic < M.size2; ic++)
                        M[ir, ic] = A[ir, ic] - B[ir, ic];
                }
                return M;
                }
        throw new ArgumentException($"Matrix -: Matrices with sizes ({A.size1}, {A.size2}) and ({B.size1}, {B.size2}) cannot be subtracted.");
}

public static matrix operator*(double a, matrix B){
	matrix M = new matrix(B.size1, B.size2);
        for(int ir=0; ir < M.size1; ir++){
        for(int ic=0; ic < M.size2; ic++)
        	M[ir, ic] = a*B[ir, ic];
        }
        return M;
}

public static matrix operator*(matrix A, double b){
	return b*A;
}

public static double[] operator*(matrix A, double[] B){
        if(A.size2 == B.Length){
                double[] C = new double[A.size1];
		for(int i=0; i < A.size1; i++){
                for(int j=0; j < A.size2; j++)
                        C[i] += A[i,j]*B[j];
                }
                return C;
                }
        throw new ArgumentException($"Matrix*: Matrix with size ({A.size1}, {A.size2}) cannot be multiplied with array of length {B.Length}.");
}

public static double[] operator*(double[] B, matrix A){
	if(A.size2 == B.Length){
                double[] C = new double[A.size2];
                for(int i=0; i < A.size2; i++){
                for(int j=0; j < A.size1; j++)
                        C[i] += A[j,i]*B[j];
                }
                return C;
                }
        throw new ArgumentException($"Matrix*: Array of size {B.Length} cannot be multiplied by matrix with size ({A.size1}, {A.size2}).");
}

public static matrix operator*(matrix A, matrix B){
	if(A.size2 == B.size1){
		matrix C = new matrix(A.size1, B.size2);
		for(int i = 0; i < C.size1; i++){
		for(int j = 0; j < C.size2; j++)
			for(int k = 0; k < A.size2; k++) 
			C[i,j] += A[i,k]*B[k,j];
		}
		return C;
	}
	throw new ArgumentException($"Matrix*: Matrices with sizes ({A.size1}, {A.size2}) and ({B.size1}, {B.size2}) cannot be multiplied");
}

public matrix copy(){
	matrix c = new matrix(size1,size2);
	for(int j=0;j<size2;j++)
		for(int i=0;i<size1;i++)
			c[i,j]=this[i,j];
	return c;
}

//String methods
public void print(string s="", string format="{0,10:g3} ", System.IO.TextWriter file=null){
	if(file==null)file = System.Console.Out;
	file.WriteLine(s);
	for(int ir=0; ir < this.size1; ir++){
	for(int ic=0; ic < this.size2; ic++)
		file.Write(format, this[ir,ic]);
		file.WriteLine();
	}
}



//Diagonals
public static matrix id(int n){
	matrix M = new matrix(n);
	for(int i = 0; i < n; i++)M[i,i]=1;
	return M;
}

public static matrix diag(vector d, int j = 0){ //constructs square matrix with a given diagonal or shifted diagonal
	matrix M = new matrix(d.size + Abs(j));
	if(j >= 0)for(int i = 0; i < d.size; i++)M[i,i+j]=d[i];
        if(j < 0)for(int i = -j; i < d.size-j; i++)M[i,i+j]=d[i+j]; 
	return M;
}

//Comparison funcs.
public static bool approx(double a, double b, double acc=1e-6, double eps=1e-6){
	if(Abs(a-b)<acc)return true;
	if(Abs(a-b)/Max(Abs(a),Abs(b)) < eps)return true;
	return false;
}

public bool approx(matrix B,double acc=1e-6, double eps=1e-6){
	if(this.size1!=B.size1)return false;
	if(this.size2!=B.size2)return false;
	for(int i=0;i<size1;i++)
		for(int j=0;j<size2;j++)
			if(!approx(this[i,j],B[i,j],acc,eps))
				return false;
	return true;
}

//Dot product of arrays
public static double dot(double[] a, double[] b){
	if( a.Length == b.Length){
		double res = 0;
		for(int i = 0; i < a.Length; i++)res+=a[i]*b[i];
		return res;
		}
	throw new ArgumentException($"Dot: Cannot dot arrays with lengths {a.Length} and {b.Length}.");
}

public static double norm(double[] a){
	return Sqrt(dot(a,a));
}

//updates
public void update(vector u, vector v, double s=1){
	for(int i=0;i<size1;i++)
	for(int j=0;j<size2;j++)
		this[i,j]+=u[i]*v[j]*s;
	}

public void sym2update(vector u, vector v, double s=1){
        for(int i=0;i<size1;i++)
        for(int j=0;j<size2;j++)
                this[i,j]+=(u[i]*v[j] + u[j]*v[i])*s;
        }
}//matrix
