package treeConstructor;

import org.apache.commons.math3.analysis.UnivariateFunction;
import seqEvolution.SubstitutionMatrix;

public class LikelihoodFunction implements UnivariateFunction {

	SubstitutionMatrix matrix;
	double[][] N;

	public LikelihoodFunction(SubstitutionMatrix matrix, double replacements[][]) {
		this.N = replacements;
		this.matrix = matrix;
	}

	// Value of derivative of log likelihood function
	@Override
	public double value(double x) throws RuntimeException {

		return matrix.getSeqDerivLogLikelihood(x, this.N);
	}
}