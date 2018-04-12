package dlModel;

import org.apache.commons.math3.analysis.UnivariateFunction;

public class TimeSampler implements UnivariateFunction {
	double s;
	double t;
	double birthRate, deathRate;
	double dy;
	double phist;
	double cdfValue;

	public TimeSampler(double birthRate, double deathRate, double t, double deathProbAtY, double phist, double nleaves,
			double cdfValue) {
		this.t = t;
		this.birthRate = birthRate;
		this.deathRate = deathRate;
		this.phist = phist;
		this.dy = deathProbAtY;
		this.s = nleaves;
		this.cdfValue = cdfValue;
	}

	// Value of cdf Function at time t
	@Override
	public double value(double to) throws RuntimeException {

		double r, p0, p1, a, dnot;
		double tminusto = this.t - to;
		double P11, SubTreeProb_N, SubTreeProb_D;

		// Computing death and subtree probability for time tminusto
		r = birthRate / deathRate;
		p0 = BirthDeathProbs.P0(birthRate, deathRate, tminusto);
		p1 = BirthDeathProbs.P1(birthRate, deathRate, tminusto);
		a = 1.0 - (r * p0 * dy);
		dnot = p0 + (p1 * dy) / a;
		SubTreeProb_N = phist * (Math.pow(r * p0, s - 1) * p1 / Math.pow(a, s + 1));

		// Computing P11 on the basis of death Prob at to
		p0 = BirthDeathProbs.P0(birthRate, deathRate, to);
		p1 = BirthDeathProbs.P1(birthRate, deathRate, to);
		a = 1.0 - (r * p0 * dnot);
		P11 = p1 / a * a;

		// Computing Subtree probability in time t
		r = birthRate / deathRate;
		p0 = BirthDeathProbs.P0(birthRate, deathRate, t);
		p1 = BirthDeathProbs.P1(birthRate, deathRate, t);
		a = 1.0 - (r * p0 * dy);
		SubTreeProb_D = phist * (Math.pow(r * p0, s - 1) * p1 / Math.pow(a, s + 1));

		if (SubTreeProb_D == 0.0) {
			System.out.println("to = " + to + "\t" + "tminusto =" + tminusto + "\t" + "P0 = " + p0 + "\t"
					+ "SubTreeProb_D = " + SubTreeProb_D + "\t" + "SubTreeProb_N = " + SubTreeProb_N);
			throw new RuntimeException();
		}

		double fv = 1.0 - ((SubTreeProb_N * P11) / SubTreeProb_D);

		return (Math.log(fv) - this.cdfValue);
	}

	// This method will test the subtree probability
	public boolean SubTreeProbIsZero() {

		// Computing Subtree probability in time t
		double r = birthRate / deathRate;
		double p0 = BirthDeathProbs.P0(birthRate, deathRate, t);
		double p1 = BirthDeathProbs.P1(birthRate, deathRate, t);
		double a = 1.0 - (r * p0 * dy);

		// Probability of the Subtree within time t
		double SubTreeProb = phist * (Math.pow(r * p0, s - 1) * p1 / Math.pow(a, s + 1));

		if (SubTreeProb == 0.0)
			return true;

		return false;
	}

	// System.out.println("f( to =" + to +" )= " + fv + "\tP11 = " + P11 +
	// "\tSubTreeProb_N = " + SubTreeProb_N);
	// return fv - this.cdfValue;
	/*
	 * ///////////////////////////////////////////////////////// // Value of cdf
	 * Function at time t public double value1(double to) throws RuntimeException {
	 * 
	 * double r,p0,p1,a,dnot; double tminusto = this.t - to; double
	 * P11,SubTreeProb_N,SubTreeProb_D;
	 * 
	 * // Computing death and subtree probability for time tminusto r =
	 * birthRate/deathRate ; p0 = BirthDeathProbs.P0(birthRate, deathRate,
	 * tminusto); p1 = BirthDeathProbs.P1(birthRate, deathRate, tminusto); a = 1.0 -
	 * ( r * p0 * dy ); dnot = p0 + (p1 * dy) / a ; SubTreeProb_N = phist * (
	 * Math.pow( r * p0, s-1) * p1 / Math.pow(a,s+1) );
	 * 
	 * // Computing P11 on the basis of death Prob at to p0 =
	 * BirthDeathProbs.P0(birthRate, deathRate, to); p1 =
	 * BirthDeathProbs.P1(birthRate, deathRate, to); a = 1.0 - ( r * p0 * dnot );
	 * P11 = p1 / a * a ;
	 * 
	 * // Computing Subtree probability in time t r = birthRate/deathRate ; p0 =
	 * BirthDeathProbs.P0(birthRate, deathRate, t); p1 =
	 * BirthDeathProbs.P1(birthRate, deathRate, t); a = 1.0 - ( r * p0 * dy );
	 * SubTreeProb_D = phist * ( Math.pow( r * p0, s-1) * p1 / Math.pow(a,s+1) );
	 * 
	 * 
	 * //double fv = 1.0 - ( (SubTreeProb_N * P11) / SubTreeProb_D ) ;
	 * 
	 * double fv = (SubTreeProb_N * P11) / SubTreeProb_D ;
	 * 
	 * //System.out.println("f( to =" + to +" )= " + fv + "\tP11 = " + P11 +
	 * "\tSubTreeProb_N = " + SubTreeProb_N);
	 * 
	 * return fv; }
	 */

}
