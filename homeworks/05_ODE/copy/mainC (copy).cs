using System;
using static System.Console;
using static System.Math;

class main{
	static int Main(){
		double t0 = 0, tL = 2.3, deltar = 1; 
		vector r0 = new vector("-0.97000436 0.24308753 0.4662036850 0.4323657300 0 0 -0.93240737 -0.86473146 0.97000436 -0.24308753 0.4662036850 0.4323657300"); /* Init. cond. Format: (x1,y1,vx1,vy1,x2,y2 ...) */
		Func<double,vector,vector> f = delegate(double t, vector r){vector w = new vector(12); /* Encoding of the differential equations */
			for(int i = 0; i < 3; i++){w[4*i] = r[4*i+2]; w[4*i+1] = r[4*i+3];} /* split the 6 equations to 12 */
			for(int i = 0; i < 3; i++) /* Each planet */
			for(int k = 0; k < 3; k++){ /* each planet different from the first*/
			if(k!=i){ deltar = Pow(Pow(r[4*k] - r[4*i],2) + Pow(r[4*k+1] - r[4*i+1],2),-1.5); /* compute the distance between planets^(-3/2) */
                        w[4*i+2] += deltar*(r[4*k] - r[4*i]); 		/* compute the force */
			w[4*i+3] += deltar*(r[4*k+1] - r[4*i+1]);}}
			return w;};
		genlist<vector> rs = new genlist<vector>();
		ODE.driver(f, t0, r0, tL, ylist: rs);
		for(int i = 0; i < rs.size; i++)WriteLine($"{rs[i][0]} {rs[i][1]} {rs[i][4]} {rs[i][5]} {rs[i][8]} {rs[i][9]}");
		return 0;
	}//Main
}//main
