package proposers;

import java.util.ArrayList;

import common.LogDouble;
import common.NormalDistribution;
import common.PRNG;
import common.Pair;
import common.RealInterval;

public class NormalProposer {

	public static Pair<Double, LogDouble[]> perturbSingleValue(double value, RealInterval interval, double proposalCV,
			PRNG prng) {

		StringBuilder logStr = new StringBuilder();

		// Initialize the forward & backward move
		LogDouble forward = new LogDouble(1.0);
		LogDouble backward = new LogDouble(1.0);

		// Initiating the normal proposal distribution
		double xOld = value;
		double stdev = Math.max(Math.abs(xOld * proposalCV), 1e-16);
		NormalDistribution pd = new NormalDistribution(xOld, Math.pow(stdev, 2));

		// Perturb the parameter by sampling values from proposal dist.
		int tries = 0;
		double x = Double.NaN;
		do {
			x = pd.sampleValue(prng);
			++tries;
			if (tries > 1000) {
				// Failed proposer, i.e., when parameters could not be perturbed within 100
				// tries
				logStr.append("Normal Proposer Error: Cannot perturb this parameter within specified interval.");
				forward = new LogDouble(1.0);
				backward = new LogDouble(0.0);
				LogDouble[] q = { forward, backward };
				Pair<Double, LogDouble[]> failProposal = new Pair<Double, LogDouble[]>(x, q);
				return failProposal;
			}
		} while (!interval.isWithin(x));

		// Obtain "forward" density.
		double a = interval.getLowerBound();
		double b = interval.getUpperBound();
		double nonTails = 1.0;
		if (!Double.isInfinite(a)) {
			nonTails -= pd.getCDF(a);
		}
		if (!Double.isInfinite(b)) {
			nonTails -= (1.0 - pd.getCDF(b));
		}
		forward.mult(new LogDouble(Math.max(pd.getPDF(x) / nonTails, 0.0)));

		// Obtain "backward" density.
		stdev = Math.max(Math.abs(x * proposalCV), 1e-16);
		pd.setMean(x);
		pd.setStandardDeviation(stdev);
		nonTails = 1.0;
		if (!Double.isInfinite(a)) {
			nonTails -= pd.getCDF(a);
		}
		if (!Double.isInfinite(b)) {
			nonTails -= (1.0 - pd.getCDF(b));
		}
		backward.mult(new LogDouble(Math.max(pd.getPDF(xOld) / nonTails, 0.0)));

		// Preparing the return value.
		LogDouble[] q = { forward, backward };
		Pair<Double, LogDouble[]> rv = new Pair<Double, LogDouble[]>(x, q);
		return rv;
	}

	// Note: this method will perturb in the passing array.
	public static LogDouble[] perturbMultiValues(double[] values, ArrayList<Integer> pIndices, RealInterval interval,
			double proposalCV, PRNG prng, boolean debugMode) {

		// Size of perturbation
		int noPerturb = pIndices.size();

		// Initialize the forward & backward move
		LogDouble forward = new LogDouble(1.0);
		LogDouble backward = new LogDouble(1.0);
		for (int i = 0; i < noPerturb; i++) {

			// Index and value to be perturbed.
			int idx = pIndices.get(i);
			double value = values[idx];

			// Initiating the normal proposal distribution
			double xOld = value;
			double stdev = Math.max(Math.abs(xOld * proposalCV), 1e-16);
			NormalDistribution pd = new NormalDistribution(xOld, Math.pow(stdev, 2));

			// Perturb the parameter by sampling values from proposal dist.
			int tries = 0;
			double x = Double.NaN;
			do {
				x = pd.sampleValue(prng);
				++tries;
				if (tries > 1000) {
					// Failed proposer, i.e., when parameters could not be perturbed within 100
					// tries
					if (debugMode)
						System.out.println(
								"Normal Proposer Error: Cannot perturb this parameter within specified interval.");

					forward = new LogDouble(1.0);
					backward = new LogDouble(0.0);
					LogDouble[] q = { forward, backward };
					return q;
				}
			} while (!interval.isWithin(x));

			// Changing old value
			// perturbedValues[idx] = x;
			values[idx] = x;

			// Obtain "forward" density.
			double a = interval.getLowerBound();
			double b = interval.getUpperBound();
			double nonTails = 1.0;
			if (!Double.isInfinite(a)) {
				nonTails -= pd.getCDF(a);
			}
			if (!Double.isInfinite(b)) {
				nonTails -= (1.0 - pd.getCDF(b));
			}
			forward.mult(new LogDouble(Math.max(pd.getPDF(x) / nonTails, 0.0)));

			// Obtain "backward" density.
			stdev = Math.max(Math.abs(x * proposalCV), 1e-16);
			pd.setMean(x);
			pd.setStandardDeviation(stdev);
			nonTails = 1.0;
			if (!Double.isInfinite(a)) {
				nonTails -= pd.getCDF(a);
			}
			if (!Double.isInfinite(b)) {
				nonTails -= (1.0 - pd.getCDF(b));
			}
			backward.mult(new LogDouble(Math.max(pd.getPDF(xOld) / nonTails, 0.0)));
		}

		// Preparing the return value.
		LogDouble[] q = { forward, backward };

		return q;
	}

}

// if(x<lbound)
// System.out.println(value);

// double lbound = interval.getLowerBound();
// if(value<lbound)
// System.out.println(value);

/*
 * double [] perturbedValues = Arrays.copyOf(values, values.length); LogDouble
 * [] q = { forward, backward }; Pair<double [],LogDouble[]> rv = new
 * Pair<double [],LogDouble[]>(perturbedValues,q); return rv; Pair<double
 * [],LogDouble[]>
 */
