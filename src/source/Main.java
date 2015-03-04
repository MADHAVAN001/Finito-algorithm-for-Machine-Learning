package source;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import Jama.Matrix;

public class Main{
	public static void main(String[]args)
	{
		int m = 10000;
		int n = 10000;
		final double mu = 10;
		String opfile = "finitoCl1-m.txt";
		final double sigma = 25;
		Matrix m_op = new Matrix(300,1);
		Matrix time = new Matrix(300,1);
		for(n = 10000;n<=10000;n+=200)
		{
			int num_cycles = 10000;

		Matrix X = new Matrix(m, n);
		Matrix Y = new Matrix(m, 1);
		/*for (int i = 0; i < X.getRowDimension(); i++)
			for (int j = 0; j < X.getColumnDimension(); j++) {
				double random = Math.random();
				double s = random - (int) random;
				X.set(i, j, mu + sigma * s);
			}
		for (int i = 0; i < Y.getRowDimension(); i++)
			for (int j = 0; j < Y.getColumnDimension(); j++) {
				double random = Math.random();
				double s = random - (int) random;
				Y.set(i, j, mu + sigma * s);
			}*/
		String Path = "Data1.txt";
		/*csvfileinput.input(X, Y, Path);
		for (int y = 0; y < X.getRowDimension(); y++)
			for (int t = 0; t < X.getColumnDimension(); t++) {
				double rand = Math.random();
				X.set(y, t, X.get(y, t) + rand * 0.01);
			}
			*/
		String Path1 = "Data/Day";
		datafileinput.input(X, Y, Path1);
		//Matrix beta_accurate = ((((X.transpose()).times(X)).inverse()).times(X
				//.transpose())).times(Y);
		Matrix beta = new Matrix(m,num_cycles);
		long startTime = System.currentTimeMillis();
		beta = finito(X,Y,num_cycles);
		long stopTime = System.currentTimeMillis();
	    long elapsedTime = stopTime - startTime;
		Matrix performance = new Matrix(num_cycles, 1);
		//Matrix abs_performance = (((X.times(beta_accurate).minus(Y))
				//.transpose()).times(X.times(beta_accurate).minus(Y)));
		performance = mse_performance.mse_problem(beta, X, Y,
				num_cycles);
		Matrix abs_performance = new Matrix(num_cycles,1);
		 //Matrix beta_accurate = new Matrix(m,1);
		 for(int tmp = 0;tmp<abs_performance.getRowDimension();tmp++)
		 {
			 abs_performance.set(tmp, 0, 0);
		 }
		/*
		 * for(int i = 0;i<beta_cd.getRowDimension();i++) { for(int j =
		 * 0;j<beta_cd.getColumnDimension();j++)
		 * System.out.printf("     "+beta_cd.get(i, j)); System.out.println(); }
		 */
		 m_op.set((n-384)/200, 0, performance.get(num_cycles-1, 0));
			time.set((n-384)/200, 0, elapsedTime);
			}
			try {
				FileWriter fw = new FileWriter(opfile);
				for (int j = 0; j < 300; j++) {

					double diff = m_op.get(j, 0)
							/ m;
					fw.append(String.valueOf(diff));
					fw.append(',');

					fw.append(String.valueOf(time.get(j, 0)));
					fw.append('\n');
				}
				fw.flush();
				fw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			}
	
	public static Matrix finito(Matrix X,Matrix Y, int num_cycles)
	{
		double alpha = Math.pow(10, -10.50);
		Matrix beta = new Matrix(X.getColumnDimension(),num_cycles);
		for(int i = 0;i<beta.getRowDimension();i++)
			for(int j = 0;j<beta.getColumnDimension();j++)
				beta.set(i, j, 0);
		Matrix g = new Matrix(X.getColumnDimension(),num_cycles);
		for(int i = 0;i<g.getRowDimension();i++)
			for(int j = 0;j<g.getColumnDimension();j++)
				g.set(i, j, 0);
		for(int i = 1;i<num_cycles;i++)
		{
			System.out.println("Cycle: " +i);
			int random;
			random = (int)Math.floor(Math.random()*X.getRowDimension());
			//System.out.println(random);
			Matrix sub_j = new Matrix(1,X.getColumnDimension());
			for(int k = 0;k<X.getColumnDimension();k++)
				sub_j.set(0, k, X.get(random,k));
			Matrix beta_cycle = new Matrix(X.getColumnDimension(),1);
			for(int k = 0;k<X.getColumnDimension();k++)
				beta_cycle.set(k, 0, beta.get(k,i-1));
			Matrix Yj = new Matrix(1,1);
			Yj.set(0,0,Y.get(random, 0));
			double tmp_i = (sub_j.times(beta_cycle).minus(Yj)).get(0, 0);
			Matrix gj = (sub_j.transpose()).times(tmp_i);
			for(int k = 0;k<gj.getRowDimension();k++)
				g.set(k, i, gj.get(k, 0));
			Matrix sum = new Matrix(g.getRowDimension(),1);
			for(int k = 0;k<g.getRowDimension();k++)
			{
				double tmp = 0;
				for(int l = 0;l<g.getColumnDimension();l++)
				{
					tmp+=g.get(k, l);
				}
				sum.set(k, 0, tmp);
			}
			Matrix tmp_beta;
			tmp_beta = beta_cycle.minus(sum.times(alpha/X.getRowDimension()));
			for(int k = 0;k<X.getColumnDimension();k++)
			{
				beta.set(k, i, tmp_beta.get(k, 0));
			}
		}
		return beta;
	}
}