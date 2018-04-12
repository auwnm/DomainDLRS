package common;

import rateModel.Normal;

/**
 * Represents a 1-D normal distribution, "Gaussian", N(m,v).
 * <p/>
 * It is possible to let the distribution rely on two floating point state
 * parameters and act as a <code>Dependent</code>.
 * 
 * @author Martin Linder.
 * @author Bengt Sennblad.
 * @author Joel Sjöstrand.
 */
public class NormalDistribution {

	/** Mean parameter. Null if not used. */
	protected Double meanParam;

	/** CV parameter. Null if not used. */
	protected Double cvParam;

	/** Mean value. Should reflect parameter's value. */
	protected double mean;

	/** Variance. Should reflect parameter's value. */
	protected double var;

	/** Standard deviation. Should reflect parameter's value. */
	protected double stdev;

	/**
	 * Constant term, given the variance, in the log density function: -0.5 * ln(2 *
	 * PI * var).
	 */
	protected double logDensFact;

	/**
	 * Constructor for when the distribution does not rely on state parameters.
	 * 
	 * @param mean
	 *            distribution mean.
	 * @param var
	 *            distribution variance.
	 */
	public NormalDistribution(double mean, double var) {
		if (var <= 0.0) {
			throw new IllegalArgumentException("Cannot have non-positive variance for normal distribution.");
		}
		this.meanParam = null;
		this.cvParam = null;
		this.mean = mean;
		this.var = var;
		this.stdev = Math.sqrt(var);
		this.update();
	}

	/**
	 * Constructor for when the distribution relies on state parameters. The
	 * distribution will add itself as a child dependent of the parameters. For MCMC
	 * mixing purposes, mean and CV has been selected as parameterisation.
	 * 
	 * @param mean
	 *            the mean parameter.
	 * @param cv
	 *            the CV parameter.
	 */
	public NormalDistribution(Double mean, Double cv) {
		this.meanParam = mean;
		this.cvParam = cv;
		this.update();
	}

	public String getName() {
		return "Normal distribution";
	}

	public int getNoOfParameters() {
		return 2;
	}

	public int getNoOfDimensions() {
		return 1;
	}

	public double getPDF(double x) {
		return Math.exp(-0.5 * Math.pow(x - this.mean, 2) / this.var + this.logDensFact);
	}

	public double getCDF(double x) {
		x = (x - this.mean) / this.stdev;
		return Normal.cdf(x);
	}

	public double getProbability(double a, double b) {
		return (this.getCDF(b) - this.getCDF(a));
	}

	public double getMean() {
		return this.mean;
	}

	public double getMedian() {
		return this.mean;
	}

	public double getStandardDeviation() {
		return this.stdev;
	}

	public double getVariance() {
		return this.var;
	}

	public double getCV() {
		return (this.stdev / Math.abs(this.mean));
	}

	/**
	 * Updates the distribution.
	 */
	public void update() {
		if (this.meanParam != null) {
			this.mean = this.meanParam;
			this.stdev = this.cvParam * this.mean;
			this.var = Math.pow(this.stdev, 2);
		}
		this.logDensFact = -0.5 * Math.log(2 * Math.PI * this.var);
	}

	@Override
	public String toString() {
		return "N(" + this.mean + ", " + this.var + ')';
	}

	public void setMean(double mean) {
		this.mean = mean;
		if (this.meanParam != null) {
			this.meanParam = mean;
		}
	}

	public void setStandardDeviation(double stdev) {
		if (stdev <= 0.0) {
			throw new IllegalArgumentException("Cannot have non-positive standard deviation for normal distribution.");
		}
		this.stdev = stdev;
		this.var = Math.pow(stdev, 2);
		if (this.cvParam != null) {
			this.cvParam = (stdev / Math.abs(this.mean));
		}
		this.logDensFact = -0.5 * Math.log(2 * Math.PI * this.var);
	}

	public void setVariance(double var) {
		if (var <= 0.0) {
			throw new IllegalArgumentException("Cannot have non-positive variance for normal distribution.");
		}
		this.var = var;
		this.stdev = Math.sqrt(var);
		if (this.cvParam != null) {
			this.cvParam = (this.stdev / Math.abs(this.mean));
		}
		this.logDensFact = -0.5 * Math.log(2 * Math.PI * this.var);
	}

	public double getQuantile(double p) {
		return (Normal.quantile(p) * this.stdev + this.mean);
	}

	public double sampleValue(PRNG prng) {
		double x = prng.nextGaussian();
		// No bounds checking for within representable range...
		return (x * this.stdev + this.mean);
	}

	public double getMode() {
		return this.mean;
	}

}
