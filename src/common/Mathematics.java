package common;

//Mainly derived from Matt Rasmussen c++ code

public class Mathematics {

	public static int choose(int n, int k) {
		return (int) (fchoose(n, k) + .5);
	}

	public static double fchoose(int n, int k) {
		if (n < 0 || k < 0 || k > n)
			return 0;

		// optimization for speed
		if (k > n / 2)
			k = n - k;

		double t = 1.0;
		double m = n;
		for (double i = 1; i <= k; i++)
			t *= (m - i + 1) / i;
		return t;
	}

	// The Java Math library function Math.random() generates a double value in the
	// range [0,1).
	// Notice this range does not include the 1.
	// In order to get a specific range of values first,
	// you need to multiply by the magnitude of the range of values you want
	// covered.

	public static int irand(int max) {
		int i = (int) (Math.random() * max);
		return (i == max) ? max - 1 : i;
	}

	public static int minIndex(double... ds) {
		int idx = -1;
		double d = Double.POSITIVE_INFINITY;
		for (int i = 0; i < ds.length; i++)
			if (ds[i] < d) {
				d = ds[i];
				idx = i;
			}
		return idx;
	}

	// This method will prepare cumulative weight distribution
	public static double[] prepareCumulative(double[] array) {
		double tot = 0.0;
		double[] cumArray = new double[array.length];
		for (int i = 0; i < array.length; ++i) {
			tot += array[i];
			cumArray[i] = tot;
		}
		for (int i = 0; i < array.length; ++i)
			cumArray[i] /= tot;

		cumArray[cumArray.length - 1] = 1.0;
		return cumArray;
	}

	/*
	 * // distributions
	 * 
	 * // probability density distribution of the Poisson double poisson(int x,
	 * float lambda) { if (x < 0 || lambda <= 0) return 0;
	 * 
	 * double a = 0; for (double i=1 ; i<x+1 ; i+=1.0) a += Math.log(lambda/i);
	 * 
	 * return Math.exp((-lambda + a) ); }
	 * 
	 * 
	 * // Normal distribution: mu is the mean, and sigma is the standard deviation
	 * public static double normalvariate(float mu, float sigma) { // Uses Kinderman
	 * and Monahan method. Reference: Kinderman, // A.J. and Monahan, J.F.,
	 * "Computer generation of random // variables using the ratio of uniform
	 * deviates", ACM Trans // Math Software, 3, (1977), pp257-260.
	 * 
	 * static double NV_MAGICCONST = (double) 4.0 * Math.exp(-0.5) / Math.sqrt(2.0);
	 * float u1, u2, z, zz;
	 * 
	 * do { u1 = frand(); u2 = 1.0 - frand(); z = NV_MAGICCONST*(u1-0.5)/u2; zz =
	 * z*z/4.0; } while (zz > -log(u2));
	 * 
	 * return mu + z*sigma; }
	 * 
	 * 
	 * } // extern "C"
	 * 
	 * 
	 * //===========================================================================
	 * == // sorting
	 * 
	 * // Invert a permutation void invertPerm(int *perm, int *inv, int size) { for
	 * (int i=0; i<size; i++) inv[perm[i]] = i; }
	 * 
	 * 
	 * 
	 * //===========================================================================
	 * == // input/output
	 * 
	 * void printIntArray(int *array, int size) { for (int i=0; i<size; i++)
	 * printf("%d ", array[i]); printf("\n"); }
	 * 
	 * void printFloatArray(float *array, int size) { for (int i=0; i<size; i++)
	 * printf("%f ", array[i]); printf("\n"); }
	 * 
	 * 
	 * 
	 * // simple random numbers
	 * //===========================================================================
	 * ==
	 * 
	 * float static frand(float max=1.0) { return rand() / float(RAND_MAX) * max; }
	 * 
	 * public static float frand(float min, float max) { return min + (rand() /
	 * float(RAND_MAX) * (max-min));
	 * 
	 * }
	 * 
	 * public static int irand(int max) { const int i = int(rand() / float(RAND_MAX)
	 * * max); return (i == max) ? max - 1 : i; }
	 * 
	 * public static int irand(int min, int max) { const int i = min + int(rand() /
	 * float(RAND_MAX) * (max - min)); return (i == max) ? max - 1 : i; }
	 * 
	 */

}
