package mcmc;

import common.LogDouble;
import common.RealInterval;

public class RealParameterUniformPrior {

	private double param; // Parameter for prior.
	public RealInterval interval; // Allowed interval.
	private boolean doUseActual; // Returns actual prior probability if true.
	private LogDouble priorProbability; // Prior.

	// Constructor.
	public RealParameterUniformPrior(RealInterval priorInterval, boolean doUseActual) {
		this.interval = priorInterval;
		this.doUseActual = doUseActual;
	}

	// This method will return the prior probability of parameter
	public LogDouble getPriorProbability(double param) {
		this.param = param;
		this.priorProbability = new LogDouble(1.0);
		if (this.doUseActual) {
			if (!this.interval.isWithin(this.param))
				this.priorProbability = new LogDouble(0.0);
			else
				this.priorProbability.mult(1.0 / this.interval.getWidth());
		} else {
			if (!this.interval.isWithin(this.param))
				this.priorProbability = new LogDouble(0.0);
		}
		return this.priorProbability;
	}

	// If doUseActual = false and parameter current value is within the interval,
	// e.g. [A,B], 1 is returned.
	// If doUseActual = true and parameter current value is within the interval,
	// 1/(B-A) is returned.
	// else 0 probability will be returned.
	public void setActualPriorProbability(boolean doUseActual) {
		this.doUseActual = doUseActual;
	}

	public String getModelName() {
		return "RealParameterUniformPrior";
	}

}
