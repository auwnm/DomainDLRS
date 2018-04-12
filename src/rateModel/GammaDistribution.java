package rateModel;

import common.LogDouble;

/**
 * Represents a 1-D gamma distribution Gamma(k,theta). It is possible to let the
 * distribution rely on two floating point state parameters and act as a
 * <code>Dependent</code>.
 * <p/>
 * Beware of the alternative parametrisation that is sometimes used. We use
 * <i>shape</i> parameter k and <i>scale</i> parameter theta. Sometimes, one
 * encounters Gamma(alpha,beta), with shape parameter alpha=k and <i>rate</i>
 * parameter beta=1/theta.
 *
 * @author Joel Sjöstrand.
 */
public class GammaDistribution {

	private double mean; // Mean state parameter. 0 if not used.
	private double cv; // CV state parameter. 0 if not used.
	protected double k; // Shape parameter. Should reflect DoubleParameters' values.
	protected double theta; // Scale parameter. Should reflect DoubleParameters' values. */
	protected double c; // For speed, holds -ln(G(k)*theta^k), where G() is the gamma function. */

	/**
	 * Constructor for when the distribution relies on state parameters. The
	 * distribution will add itself as a child dependent of the parameters. For MCMC
	 * mixing purposes, mean and CV has been selected as parameterisation.
	 * 
	 * @param mean
	 *            the mean.
	 * @param cv
	 *            the CV.
	 */
	public GammaDistribution(double mean, double cv) {
		if (mean <= 0.0 || cv <= 0.0) {
			throw new IllegalArgumentException("Cannot have non-positive mean or cv for gamma distribution.");
		}
		this.k = 0;
		this.theta = 0;
		this.mean = mean;
		this.cv = cv;
		this.update();
	}

	/*
	 * Constructor for when the distribution does not rely on state parameters. Note
	 * the parametrisation where <i>scale</i> parameter theta=1/beta, for
	 * <i>rate</i> parameter beta.
	 * 
	 * @param k distribution shape parameter.
	 * 
	 * @param theta distribution scale parameter.
	 */
	public GammaDistribution(Double k, Double theta) {
		if (k <= 0.0 || theta <= 0.0) {
			throw new IllegalArgumentException("Cannot have non-positive shape or scale for gamma distribution.");
		}
		this.mean = 0;
		this.cv = 0;
		this.k = k.doubleValue();
		this.theta = theta.doubleValue();
		this.update();
	}

	/**
	 * Updates the distribution following, e.g., a change of the underlying
	 * DoubleParameters.
	 */
	private void update() {
		if (this.mean != 0 && this.cv != 0) {
			double cv2 = this.cv * this.cv;
			this.k = 1.0 / cv2;
			this.theta = this.mean * cv2;
		}
		this.c = -this.k * Math.log(this.theta) - Gamma.lnGamma(this.k);
	}

	public double getPDF(double x) {
		return Math.exp((this.k - 1.0) * Math.log(x) - x / this.theta + this.c);
	}

	/**
	 * For speed reasons, provides access to the probability density function f(x)
	 * for a specified value x, returned as a <code>Probability</code>.
	 * 
	 * @param x
	 *            the value where to evaluate density.
	 * @return the probability density, f(x).
	 */
	public LogDouble getPDFAsProbability(double x) {
		return new LogDouble((this.k - 1.0) * Math.log(x) - x / this.theta + this.c, 1);
	}

	public double getCDF(double x) {
		return Gamma.gammaCDF(x, this.k, this.theta);
	}

	public double getQuantile(double p) {
		return Gamma.gammaQuantile(p, this.k, this.theta);
	}

	public double getProbability(double a, double b) {
		return (this.getCDF(b) - this.getCDF(a));
	}

	public double getMean() {
		return (this.k * this.theta);
	}

	/**
	 * Sets the mean by changing the shape parameter but not scale parameter.
	 * 
	 * @param mean
	 *            the new mean.
	 */

	public void setMean(double mean) {
		if (mean <= 0.0) {
			throw new IllegalArgumentException("Cannot set non-positive mean on gamma distribution.");
		}
		this.k = mean / this.theta;
		this.c = -this.k * Math.log(this.theta) - Gamma.lnGamma(this.k);
		if (this.mean != 0) {
			synchMeanAndCV();
		}
	}

	public double getMedian() {
		return this.getQuantile(0.5);
	}

	public double getStandardDeviation() {
		return Math.sqrt(this.k) * this.theta;
	}

	/**
	 * Sets the standard deviation by changing the shape parameter but not scale
	 * parameter.
	 * 
	 * @param stdev
	 *            the new standard deviation.
	 */

	public void setStandardDeviation(double stdev) {
		if (stdev <= 0.0) {
			throw new IllegalArgumentException("Cannot set non-positive variance on gamma distribution.");
		}
		this.k = Math.pow(stdev / this.theta, 2);
		this.c = -this.k * Math.log(this.theta) - Gamma.lnGamma(this.k);
		if (this.mean != 0) {
			synchMeanAndCV();
		}
	}

	public double getVariance() {
		return (this.k * this.theta * this.theta);
	}

	/**
	 * Sets the variance by changing the shape parameter but not scale parameter.
	 * 
	 * @param var
	 *            the new variance.
	 */

	public void setVariance(double var) {
		if (var <= 0.0) {
			throw new IllegalArgumentException("Cannot set non-positive variance on gamma distribution.");
		}
		this.k = var / (this.theta * this.theta);
		this.c = -this.k * Math.log(this.theta) - Gamma.lnGamma(this.k);
		if (this.mean != 0) {
			synchMeanAndCV();
		}
	}

	/**
	 * Adjusts DoubleParameters' values in accordance with the current internal
	 * values following a change of the shape parameter.
	 */
	private void synchMeanAndCV() {
		this.mean = this.k * this.theta;
		this.cv = 1.0 / Math.sqrt(this.k);
	}

	public double getCV() {
		return (1.0 / Math.sqrt(k));
	}

	public double sampleValue(double randno) {
		return Gamma.gammaQuantile(randno, this.k, this.theta);
	}

	@Override
	public String toString() {
		return "Gamma(" + this.k + ", " + this.theta + ')';
	}

	public double getMode() {
		if (this.k >= 1.0) {
			return ((this.k - 1.0) * this.theta);
		}
		throw new UnsupportedOperationException("Cannot compute gamma distribution mode when shape parameter k < 1.");
	}

	/**
	 * Returns the shape parameter.
	 * 
	 * @return the shape parameter k, a.k.a alpha.
	 */
	public double getShape() {
		return this.k;
	}

	/**
	 * Returns the scale parameter theta.
	 * 
	 * @return the scale parameter theta, a.k.a beta^-1.
	 */
	public double getScale() {
		return this.theta;
	}

}
