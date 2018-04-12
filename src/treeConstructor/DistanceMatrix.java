package treeConstructor;

import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.analysis.solvers.UnivariateSolver;

import seqEvolution.MSAData;
import seqEvolution.SequenceType;
import seqEvolution.SubstitutionMatrix;

public class DistanceMatrix {

	public static final int PID = 0;
	public static final int JC69 = 1;
	public static final int JCAA = 2;
	public static final int KIMURA = 3;
	public static final int POISSON = 4;
	public static final int MLE = 5;

	public static double[][] compute(MSAData msa, SubstitutionMatrix sm, int type) throws RuntimeException {
		if (type == PID)
			return computeID(msa);
		else if (type == JC69)
			return computeJC69(msa);
		else if (type == JCAA)
			return computeTajimaNei(msa);
		else if (type == KIMURA)
			return computeKIMURA(msa);
		else if (type == POISSON)
			return computePoisson(msa);
		else if (type == MLE)
			return computeMLE_BR(msa, sm);
		else {
			System.out.println("Error: Specified model is not available.");
			return null;
		}
	}

	// Distance matrix based on Maximum Likelihood Estimate using Brent Method.
	// In numerical analysis, Brent's method is a complicated but popular
	// root-finding algorithm combining
	// the bisection method, the secant method and inverse quadratic interpolation.
	public static double[][] computeMLE_BR(MSAData msa, SubstitutionMatrix sm) throws RuntimeException {

		// Initializing it with Jukes Cantor Jukes-Cantor correction distance
		double[][] distance = computeTajimaNei(msa);

		double mleDist = 0;
		final double minValue = 0.00001;
		final double maxValue = SubstitutionMatrix.MAX_MARKOV_TIME;
		final int maxIterations = 1000;
		final double relativeAccuracy = 1.0e-12;
		final double absoluteAccuracy = 1.0e-8;

		int m = msa.getNoOfSequences();

		for (int i = 0; i < (m - 1); i++) {
			for (int j = i; j < m; j++) {
				if (j != i) {
					double N[][] = countReplacements(msa.dataAsStrings[i], msa.dataAsStrings[j], msa.getSequenceType());
					LikelihoodFunction lFunction = new LikelihoodFunction(sm, N);
					if ((lFunction.value(minValue) < 0 && lFunction.value(maxValue) < 0)
							|| (lFunction.value(minValue) > 0 && lFunction.value(maxValue) > 0)) {
						System.out.println("Warning: Root is not bracketd within specified interval for " + (i + 1)
								+ " and " + (j + 1));
					} else {
						UnivariateSolver nonBracketing = new BrentSolver(relativeAccuracy, absoluteAccuracy);
						mleDist = nonBracketing.solve(maxIterations, lFunction, minValue, maxValue);
						distance[i][j] = mleDist;
						distance[j][i] = distance[i][j];
					}
				}
			}
		}
		return distance;
	}

	// Distance matrix based on Maximum Likelihood Estimate using Brent Method.
	// In numerical analysis, Brent's method is a complicated but popular
	// root-finding algorithm combining
	// the bisection method, the secant method and inverse quadratic interpolation.
	public static double computeMLE_BR(String seq1, String seq2, SubstitutionMatrix sm) throws RuntimeException {

		// Initializing it with Jukes Cantor Jukes-Cantor correction distance
		double distance = Math.random();

		double mleDist = 0;
		final double minValue = 1.0e-12;
		final double maxValue = SubstitutionMatrix.MAX_MARKOV_TIME;
		final int maxIterations = 1000;
		final double relativeAccuracy = 1.0e-12;
		final double absoluteAccuracy = 1.0e-8;

		double N[][] = countReplacements(seq1, seq2, sm.getSequenceType());

		LikelihoodFunction lFunction = new LikelihoodFunction(sm, N);

		double minFun = lFunction.value(minValue);
		double maxFun = lFunction.value(maxValue);

		// if( (lFunction.value(minValue)<0 && lFunction.value(maxValue)<0) ||
		// (lFunction.value(minValue)>0 && lFunction.value(maxValue)>0) ) {
		if (minFun * maxFun > 0) {
			System.out.println("Warning: Root is not bracketd within specified interval.");
		} else {
			UnivariateSolver nonBracketing = new BrentSolver(relativeAccuracy, absoluteAccuracy);
			mleDist = nonBracketing.solve(maxIterations, lFunction, minValue, maxValue);
			distance = mleDist;
		}

		return distance;
	}

	// Derived from Fast Phylo code
	// Available at
	// http://sourceforge.net/p/fastphylo/code/HEAD/tree/trunk/src/c++/programs/fastprot/MaximumLikelihood.cpp#l7
	// Computes the distance with the maximum likelihood using Newton-Rhapson
	public static double[][] computeMLE_NR(MSAData msa, SubstitutionMatrix sm) {

		// Initializing it with kimura distance
		double[][] distance = computeKIMURA(msa);
		int m = msa.getNoOfSequences();

		for (int i = 0; i < (m - 1); i++) {
			for (int j = i; j < m; j++) {
				if (j != i) {

					double N[][] = countReplacements(msa.dataAsStrings[i], msa.dataAsStrings[j], msa.getSequenceType());
					LikelihoodFunction lFunction = new LikelihoodFunction(sm, N);

					double t = distance[i][j];
					if (t == 0)
						t = 1;

					// Newton-Rhapson
					double delta, tol;
					delta = tol = 0.001;
					int maxIterations = 50;
					double l_d = lFunction.value(t);

					for (int k = 0; k < maxIterations; k++) {

						if (Math.abs(l_d) < tol) { // If the derivative is small enough
							break;
						}
						double l_new = lFunction.value(t + delta);
						double deriv = (l_new - l_d) / delta;
						t = t - l_d / deriv;
						if (t < 1) { // 1 is the smallest possible distance
							t = 1;
							break;
						}
						if (t > 500) { // 500 is infinity
							t = 500;
							break;
						}
						if (Math.abs(l_d) < Math.abs(l_new)) { // The derivative is getting larger..
							break;
						}
						l_d = l_new;
					}
					distance[i][j] = t;
					distance[j][i] = distance[i][j];
				}
			}
		}
		return distance;
	}

	// Computes the identity based distance
	public static double[][] computeID(MSAData msa) {

		int m = msa.getNoOfSequences();
		double[][] distance = new double[m][m];

		for (int i = 0; i < (m - 1); i++) {
			for (int j = i; j < m; j++) {
				if (j == i) {
					distance[i][i] = 0;
				} else {

					distance[i][j] = 1.0
							- countIdentities(msa.dataAsStrings[i], msa.dataAsStrings[j], msa.getSequenceType());
					distance[j][i] = distance[i][j];
				}
			}
		}

		return distance;
	}

	// Poisson correction
	// The distance between two amino acid sequences is computed starting from the
	// assumption that
	// the rate of amino acid substitution at each site follows the poisson
	// distribution (e.g. Zuckerkandl and Pauling, 1965; Dickerson, 1971):
	// d = -ln ( 1 - p )
	// Where p is the fraction of different amino acids between two sequences
	// (dissimilarity).
	// Ref:
	// http://bioinformatics.psb.ugent.be/downloads/psb/Userman/treecon_distance.html
	public static double[][] computePoisson(MSAData msa) {

		int m = msa.getNoOfSequences();
		double[][] distance = new double[m][m];

		for (int i = 0; i < (m - 1); i++) {
			for (int j = i; j < m; j++) {
				if (j == i) {
					distance[i][i] = 0;
				} else {
					double p = 1.0 - countIdentities(msa.dataAsStrings[i], msa.dataAsStrings[j], msa.getSequenceType());
					distance[i][j] = -1.0 * Math.log(1.0 - p);
					distance[j][i] = distance[i][j];
				}
			}
		}
		return distance;

	}

	// Evolutionary distances of Jukes and Cantor (1969)
	// d = -(3/4) * ln ( 1 - (4/3) * p )
	// Where p is the fraction of different amino acids between two sequences
	// (dissimilarity).
	// Ref:
	// http://bioinformatics.psb.ugent.be/downloads/psb/Userman/treecon_distance.html
	public static double[][] computeJC69(MSAData msa) {

		int m = msa.getNoOfSequences();
		double[][] distance = new double[m][m];
		double b = 3.0 / 4.0;

		for (int i = 0; i < (m - 1); i++) {
			for (int j = i; j < m; j++) {
				if (j == i) {
					distance[i][i] = 0;
				} else {
					double p = 1.0 - countIdentities(msa.dataAsStrings[i], msa.dataAsStrings[j], msa.getSequenceType());
					if (p >= b)
						p = b - 0.00001;
					distance[i][j] = -1.0 * b * Math.log(1.0 - (1.0 / b) * p);
					distance[j][i] = distance[i][j];
				}
			}
		}
		return distance;

	}

	// (Jukes and Cantor 1969) approximation for amino acid sequences
	// d = -b * ln ( 1 - (1/b) * p )
	// Where p is the fraction of different amino acids between two sequences
	// (dissimilarity).
	// and b= 0.95 (i.e. 19/20) for amino acid replacements (Tajima and Nei, 1984).
	// Ref:
	// http://bioinformatics.psb.ugent.be/downloads/psb/Userman/treecon_distance.html
	public static double[][] computeTajimaNei(MSAData msa) {

		int m = msa.getNoOfSequences();
		double[][] distance = new double[m][m];
		double b = 19.0 / 20.0;

		for (int i = 0; i < (m - 1); i++) {
			for (int j = i; j < m; j++) {
				if (j == i) {
					distance[i][i] = 0;
				} else {
					double p = 1.0 - countIdentities(msa.dataAsStrings[i], msa.dataAsStrings[j], msa.getSequenceType());
					distance[i][j] = -1.0 * b * Math.log(1.0 - (1.0 / b) * p);
					distance[j][i] = distance[i][j];
				}
			}
		}
		return distance;

	}

	// Calculates the Kimura corrected distance between two sequences
	// d = - ln ( 1 - p - 0.2 p^2 ).
	// Ref: http://evolution.genetics.washington.edu/phylip/doc/protdist.html
	public static double[][] computeKIMURA(MSAData msa) {

		int m = msa.getNoOfSequences();
		double[][] distance = new double[m][m];

		for (int i = 0; i < (m - 1); i++) {
			for (int j = i; j < m; j++) {
				if (j == i) {
					distance[i][i] = 0;
				} else {
					double p = 1.0 - countIdentities(msa.dataAsStrings[i], msa.dataAsStrings[j], msa.getSequenceType());
					double adj_distance = p + 0.2 * p * p;
					if (adj_distance > 0.854)
						adj_distance = 0.854;

					distance[i][j] = -1.0 * Math.log(1 - adj_distance);
					distance[j][i] = distance[i][j];
				}
			}
		}

		return distance;
	}

	// Counts all replacements from one amino acid to another in two sequences
	public static double[][] countReplacements(String s1, String s2, SequenceType st) {

		// precondition
		assert (s1.length() == s2.length());

		int alignLen = s1.length();
		double[][] N = new double[st.getAlphabetSize()][st.getAlphabetSize()];

		for (int i = 0; i < alignLen; i++) {
			int c1 = st.char2int(s1.charAt(i));
			int c2 = st.char2int(s2.charAt(i));
			if (c1 < st.getAlphabetSize() && c2 < st.getAlphabetSize())
				N[c1][c2]++;
		}
		return N;
	}

	// Counts the percentage identity, the number of amino acids that are
	// the same between two sequences, divided with the length of the sequences
	// will return The percentage identity
	public static double countIdentities(String s1, String s2, SequenceType st) {

		// precondition
		assert (s1.length() == s2.length());

		int alignLen = s1.length();
		double identicalPos = 0;

		for (int i = 0; i < alignLen; i++) {
			int c1 = st.char2int(s1.charAt(i));
			int c2 = st.char2int(s2.charAt(i));

			if (c1 < st.getAlphabetSize() && c2 < st.getAlphabetSize())
				if (c1 == c2)
					identicalPos++;
		}

		return identicalPos / alignLen;
	}

	// Simple utility method for observing the distance matrix
	public static void show(double[][] distances) {
		int rows = distances.length;
		int cols = distances[0].length;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				System.out.print(String.format("%.2f", distances[i][j]) + "\t");
			}
			System.out.println();
		}
	}

}
