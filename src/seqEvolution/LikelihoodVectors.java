package seqEvolution;

import org.ejml.data.DenseMatrix64F;

public class LikelihoodVectors {

	public DenseMatrix64F[] likelihoods;

	// Constructor
	public LikelihoodVectors(int noOfPatterns, int alphabetSize) {
		likelihoods = new DenseMatrix64F[noOfPatterns];
		for (int i = 0; i < noOfPatterns; ++i) {
			this.likelihoods[i] = new DenseMatrix64F(alphabetSize, 1);
		}
	}

}
