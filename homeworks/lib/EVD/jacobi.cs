using static System.Math;

public static class jacobi{
	public static void timesJ(matrix A, int p, int q, double theta){
		double c = Cos(theta), s = Sin(theta);
		for(int i = 0; i < A.size1; i++){
			double Aip = A[i,p], Aiq = A[i,q];
			A[i,p] = c*Aip - s*Aiq;
			A[i,q] = c*Aiq + s*Aip;
			}
	}//timesJ
	
	public static void Jtimes(matrix A, int p, int q, double theta){
                double c = Cos(theta), s = Sin(theta);
                for(int i = 0; i < A.size2; i++){
			double Api = A[p,i], Aqi = A[q,i];
                        A[p,i] = c*Api + s*Aqi;
                        A[q,i] = c*Aqi - s*Api;
                        }
	}//Jtimes
	 
	public static (vector, matrix) cyclic(matrix A){
		matrix D = A.copy(), V = matrix.id(A.size1);	
		vector w = new vector(D.size1);
		bool changed;
		do{
			changed = false; //assume succesful diagonalization until proven otherwise
			for(int p = 0; p < D.size1-1; p++) //loop over upper triangle part of matrix to do rotation
			for(int q = p+1; q < D.size2; q++){
                                double theta = 0.5*Atan2(2*D[p,q], D[q,q] - D[p,p]); //find theta
				double c = Cos(theta), s = Sin(theta);
				double new_Dpp = c*c*D[p,p] - 2*s*c*D[p,q] + s*s*D[q,q]; //find new diagonal elements
				double new_Dqq = c*c*D[q,q] + 2*s*c*D[p,q] + s*s*D[p,p];
				if(new_Dpp != D[p,p] || new_Dqq != D[q,q]){
                                	timesJ(D, p, q, theta);
                                	Jtimes(D, p ,q, -theta); //D <- J^T*D*J
                                	timesJ(V, p, q, theta); //V <- VJ
					changed = true;} //proven otherwise
				}
		}while(changed);
		for(int i = 0; i < w.size; i++)w[i] = D[i,i]; //collect the eigenvalues as the diagonal elements of D
		return (w, V);
	}//cyclic
	
	public static (vector, matrix) cyclic_opt(matrix A){
		matrix D = A.copy(), V = matrix.id(A.size1);
                vector w = new vector(D.size1);
                bool changed;
                do{
                        changed = false; //assume succesful diagonalization until proven otherwise
                        for(int p = 0; p < D.size1-1; p++) //loop over upper triangle part of matrix to do rotation
                        for(int q = p+1; q < D.size2; q++){
                                double theta = 0.5*Atan2(2*D[p,q], D[q,q] - D[p,p]); //find theta
                                double c = Cos(theta), s = Sin(theta);
                                double new_Dpp = c*c*D[p,p] - 2*s*c*D[p,q] + s*s*D[q,q]; //find new diagonal elements
                                double new_Dqq = c*c*D[q,q] + 2*s*c*D[p,q] + s*s*D[p,p];
                                if(new_Dpp != D[p,p] || new_Dqq != D[q,q]){
					changed = true;
					timesJ(V, p, q, theta); 			//V <- VJ
					D[p,q] = s*c*(D[p,p] - D[q,q]) + (c*c - s*s)*D[p,q]; //update cross point
					D[p,p] = new_Dpp; D[q,q] = new_Dqq; 		//update diagonal
					for(int i = 0; i < D.size2; i++){
					if(i!=q && i!=p){
						int Mqi = Max(q,i), mqi = Min(q,i), Mpi = Max(p,i), mpi = Min(p,i); 
						double Dpi = D[mpi, Mpi], Dqi = D[mqi,Mqi]; //set old data 
						D[mpi,Mpi] = c*Dpi - s*Dqi;			//update data
						D[mqi,Mqi] = s*Dpi + c*Dqi;}}
					}
			}
		}while(changed);
		for(int i = 0; i < w.size; i++)w[i] = D[i,i]; //collect the eigenvalues as the diagonal elements of D
                return (w, V);
	}//cyclic_opt
}//jacobi
