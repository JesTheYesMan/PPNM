using System;
using static System.Math;

namespace min{

	public class newton{
		public static readonly double epsilon = Pow(2,-26);	
        public static readonly double epsilon_c = Pow(2,-20);
		public static readonly double epsilon_4 = Pow(2,-12);
		public static readonly double lambda_min=Pow(2,-13);

        public readonly Func<vector,double> F;
        public readonly int n;
        public vector x;
        public double f;
        public int steps;
        public int f_eval;
        public bool status;

        public newton(Func<vector,double> func, vector start, double acc = 0.001, int max_steps = 9999){
            F=func; n = start.size; steps = 0; f_eval = 0;
            x = start.copy();
        	double lambda = 1, fx = 0, fy = 0;
		
			do{
				steps++;
				fx = F(x); f_eval++;
				(vector gradx, matrix H) = grad_hess(x,fx);
				if(gradx.norm() < acc) break;
				vector Dx = QRGS.solve(H, -gradx);
				lambda = 1;
				do{
					fy = F(x+lambda*Dx); f_eval++;
					if(fy < fx ) break;
					if(lambda < lambda_min) break;
					lambda /= 2;
				}while(true);
				fx = fy; x+=lambda*Dx;
			}while(steps < max_steps);
            f = fx;
            status = steps < max_steps;
        }

		(vector,matrix) grad_hess(vector x, double fx = Double.NaN){
			vector gradx = gradient(x,fx);
			return (gradx,hessian(x,fx,gradx));
		}

		vector gradient(vector x,double fx = Double.NaN){
			vector grad = new vector(n);
			if(Double.IsNaN(fx)){fx = F(x); f_eval++;}
			for(int i=0;i<n;i++){
				double dx=Abs(x[i])*epsilon;
				x[i]+=dx;
				grad[i]=(F(x)-fx)/dx; 
				f_eval++;
				x[i]-=dx;
			}
			return grad;
		}

		matrix hessian(vector x,double fx = Double.NaN,vector gradx = null){
			matrix H=new matrix(n);
            if(Double.IsNaN(fx)){fx = F(x); f_eval++;}
			if(gradx==null) gradx=gradient(x,fx);
			for(int j=0;j<n;j++){
				double dx=Abs(x[j])*epsilon_4;
				x[j]+=dx;
				vector dgrad=gradient(x)-gradx;
				for(int i=0;i<n;i++) H[i,j]=dgrad[i]/dx;
				x[j]-=dx;
			}
			return 0.5*(H+H.transpose());
		}
	}
	
	public class qnewton{ 
		static readonly double lambda_min=Pow(2,-26);
		
		public readonly Func<vector,double> F;
		public readonly int n;
		public vector x;
		public double f;
		public int steps;
		public int f_eval;
		public bool status;

		public qnewton(Func<vector,double> func, vector start, double acc = 1e-3, int max_steps = 9999){
			F=func; n = start.size; steps = 0; f_eval = 0;
			matrix B = matrix.id(n);
			vector gradf = new vector(n), Dx = new vector(n); 
			x = start.copy();
			double lambda = 1, fx = 0, fy = 0, gamma = 0, sy = 0;
			vector y = new vector(n), u = new vector(n), a = new vector(n);

			do{
				steps++;
				fx = F(x); f_eval++;
				gradf = grad(x,fx);
				Dx = -B*gradf;
				lambda = 1;
				while(true){
					fy = F(x+lambda*Dx); f_eval++;
					if(fy<fx){
						x+=lambda*Dx;
						y = grad(x,fy)-gradf;
						sy = lambda*Dx.dot(y);
						if(Abs(sy) > Pow(10,-6)){
							u = lambda*Dx - B*y;
							gamma = u.dot(y)/(2*sy);
							a = (u-gamma*lambda*Dx)/sy;
							B.sym2update(a,Dx,lambda);
						}
						break;
					}
					lambda/=2;
					if(lambda < lambda_min){
						x+=lambda*Dx;
						B = matrix.id(n);
						break;
					}
				}
			}while(gradf.norm()>acc && steps < max_steps);
			f = fx;
			status = steps < max_steps;
		}

        vector grad(vector x, double f0 = double.NaN){
        	if(f0 == double.NaN) f0 = F(x);
			double eps = Pow(2,-26);
           	vector res = new vector(n), y = x.copy();
           	for(int i=0;i<n;i++){
                double dx = Abs(x[i])*eps;
                y[i]+=dx;
                res[i] = (F(y)-f0)/dx;
				f_eval++;
                y[i]=x[i];
			}
            return res;
        }
	}
	
	public class downhill_sim{
		public readonly Func<vector,double> F;
        public readonly int n;
        public vector x;
        public double f;
        public int steps;
        public int f_eval;
        public bool status;

		public downhill_sim(Func<vector,double> func, vector start, double acc = 1e-10, int max_steps = 9999, double d = 5){
			F = func;
			steps = 0; status = false;
			simplex simp = new simplex(F,start,d:d);
			do{
				simp.update_op_vals();
				if(simp.ref_val < simp.minval){
					if(simp.exp_val < simp.ref_val)simp.expansion();
					else simp.reflection();
				}
				else{
					if(simp.ref_val < simp.maxval)simp.reflection();
					else if(simp.con_val < simp.maxval)simp.contraction();
					else simp.reduction();
				}
				steps++;
			}while(simp.std()>acc && steps < max_steps);
			status = steps < max_steps;
			x = simp.min();
			f = simp.minval;
			f_eval = simp.simp_f_eval;
		}

		class simplex{
			public readonly int dim;
			public int nmax, nmin, simp_f_eval;
			public vector[] points;
			public vector phi_vals;
			public readonly Func<vector,double> phi;
			public vector centroid;
			public double maxval, minval, ref_val, exp_val, con_val;

			public simplex(Func<vector,double> f, vector x, double d=5){
				dim = x.size; phi = f; nmax = 0; nmin = 0;
				vector[] ps = new vector[dim+1]; ps[0] = x.copy();
				for(int i=1;i<dim+1;i++){ps[i]=x.copy();ps[i][i-1]+=d;}
				points = ps; simp_f_eval = 0;
        	    phi_vals = new double[dim+1]; 
				update_phi_vals(); update_cent();
			}

			void update_phi_vals(){
				for(int i=0;i<dim+1;i++)phi_vals[i] = phi(points[i]);
				simp_f_eval += dim+1;
				nmin = 0; nmax = 0; maxval = phi_vals[0]; minval = phi_vals[0];
				for(int i=1;i<dim+1;i++){
					if(phi_vals[i]<minval){nmin=i; minval = phi_vals[i];}
					if(phi_vals[i]>maxval){nmax=i; maxval = phi_vals[i];}
				}
			}

			void update_max(){
				nmax = 0; maxval = phi_vals[0];
				for(int i=1;i<dim+1;i++){
					if(phi_vals[i]>maxval){nmax=i; maxval = phi_vals[i];}
				}
			}

			void update_cent(){
				centroid = new vector(dim);
				for(int i=0;i<dim+1;i++){
					if(i!=nmax)centroid+=points[i]/dim;
				}
			}
			
			public void update_op_vals(){
				ref_val = phi(2*centroid - points[nmax]);
				exp_val = phi(3*centroid - 2*points[nmax]);
				con_val = phi(centroid/2 - points[nmax]/2);
				simp_f_eval += 3;
			}

			public vector min(){return points[nmin];}
			
			public void reflection(){
				points[nmax] = 2*centroid - points[nmax];
                phi_vals[nmax] = ref_val;
				if(ref_val < minval){nmin=nmax; minval = ref_val;}
				update_max(); update_cent();
			}

			public void expansion(){
				points[nmax] = 3*centroid - 2*points[nmax];
                phi_vals[nmax] = exp_val;
                if(exp_val < minval){nmin=nmax; minval = exp_val;}
                update_max(); update_cent();
			}

			public void contraction(){
	            points[nmax] = centroid/2 - points[nmax]/2;
                phi_vals[nmax] = con_val;
                if(con_val < minval){nmin=nmax; minval = con_val;}
                update_max(); update_cent();
	        }
				
			public void reduction(){
				for(int i=0;i<dim+1;i++){
					if(i!=nmin)points[i]=(points[i]+points[nmin])/2;
				}
				update_phi_vals(); update_cent();
			}

			public double std(){
				double mean = 0, mean_sq = 0;
				for(int i=0;i<dim+1;i++){
					mean += phi_vals[i]/(dim+1);
					mean_sq += Pow(phi_vals[i],2)/(dim+1);
				}
				return mean_sq - mean*mean;
			}
		}
	}
	public class glo_sto{

	public readonly int dim, f_eval;
	public readonly vector x;
	public readonly double f;
	public readonly bool status;

	//constructor
	public glo_sto(Func<vector,double> φ, vector a, vector b, int seconds,double acc =1e-10,double d = 1){
		dim = a.size;
		f_eval = 0;
		f = Double.PositiveInfinity;
		x = new vector(dim);
		vector x_temp = new vector(dim);
		var start_time = DateTime.Now;
		do{
			f_eval++;
			halton(f_eval, dim, x_temp);
			for(int i=0;i<dim;i++) x_temp[i] = a[i] + x_temp[i]*(b[i]-a[i]);
			double f_temp = φ(x_temp);
			if(f_temp < f){f = f_temp; x = x_temp;}
		}while((DateTime.Now-start_time).Seconds < seconds);
		min.downhill_sim mini = new min.downhill_sim(φ,x,acc,d:d);
		x = mini.x;
		f = mini.f;
		f_eval += mini.f_eval;
		status = mini.status;
	}

	static double corput(int n, int b){
		double q=0, bk=(double)1.0/b;
		while(n>0){q+=(n%b)*bk; n/=b; bk/=b;}
		return q;
	}//corput

	static void halton(int n, int d, vector x){
		int[] basis = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293};
		if(d > basis.Length) throw new ArgumentException($"Halton: dimension too large: {d} > {basis.Length}.");
		for(int i = 0; i < d; i++) x[i] = corput(n, basis[i]);
	}//halton

	}
}

