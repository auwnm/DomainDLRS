package proposers;

import java.util.Arrays;

import Main.Parameters;

import common.LogDouble;
import common.PRNG;
import common.Pair;
import common.RealInterval;

public class RateProposer implements Proposer {

	public String name;
	public Parameters pi;
	public int type;
	public int level;
	public int domIndex;

	public PRNG prng;
	public double proposalCV;
	public RealInterval interval;
	public ProposerStatistics stats;

	public double[] ratesCache; // Rate parameter cache
	public int perturbIdx; // Index of the value going to be perturbed

	private LogDouble forwardDensity; // Forward probability density.
	private LogDouble backwardDensity; // Backward probability density.
	private boolean isEnabled; // On/off switch. By default it will be true.

	public RateProposer(int type, Parameters pi, int level, int domIndex, double proposalCV, PRNG prng) {
		this.pi = pi;
		this.level = level;
		this.type = type;
		this.domIndex = domIndex;

		if (this.level == Parameters.GENE_LEVEL)
			this.name = "geneBDRates";
		else
			this.name = (this.type == Parameters.BD_RATES) ? "domBDRatesProp_" + (this.domIndex + 1)
					: "domEdgeRatesProp_" + (this.domIndex + 1);

		this.isEnabled = true;
		this.ratesCache = null;
		this.prng = prng;
		this.proposalCV = proposalCV;
		this.interval = new RealInterval(0, Double.POSITIVE_INFINITY, true, true);
		this.stats = new ProposerStatistics(this.name);
	}

	// This method will perturb single rate parameter w.r.t. its randomly choosen
	// index
	@Override
	public boolean cacheAndPerturb() {

		// Specify the rate parameter
		double[] rates;
		if (this.level == Parameters.GENE_LEVEL)
			rates = pi.geneBDRates;
		else
			rates = (this.type == Parameters.BD_RATES) ? pi.domBDRates[this.domIndex] : pi.domEdgeRates[this.domIndex];

		// Cache the current state of rate parameter(s)
		this.ratesCache = Arrays.copyOf(rates, rates.length);

		// Randomly choose the rate parameter for perturbation
		this.perturbIdx = this.prng.nextInt(rates.length);
		Pair<Double, LogDouble[]> perturbedValue = NormalProposer.perturbSingleValue(rates[this.perturbIdx],
				this.interval, this.proposalCV, this.prng);

		// Changing the state of parameter
		rates[this.perturbIdx] = perturbedValue.first;
		this.forwardDensity = perturbedValue.second[0];
		this.backwardDensity = perturbedValue.second[1];

		return true;
	}

	// Clears the cached (previous) state of perturbed parameters when a proposed
	// state has been accepted.
	@Override
	public void clearCache() {
		if (this.stats != null) {
			String category;
			if (this.level == Parameters.GENE_LEVEL)
				category = this.perturbIdx == 0 ? "birthRate" : "deathRate";
			else {
				if (this.type == Parameters.BD_RATES)
					category = this.perturbIdx == 0 ? "birthRate" : "deathRate";
				else
					category = this.perturbIdx == 0 ? "edgeMean" : "edgeCV";
			}
			this.stats.increment(true, category);
		}
		this.ratesCache = null;
		this.perturbIdx = -1;
	}

	// Restores the cached (previous) state of perturbed parameters when a proposed
	// state has been rejected.
	@Override
	public void restoreCache() {
		if (this.stats != null) {
			String category;
			if (this.level == Parameters.GENE_LEVEL)
				category = this.perturbIdx == 0 ? "birthRate" : "deathRate";
			else {
				if (this.type == Parameters.BD_RATES)
					category = this.perturbIdx == 0 ? "birthRate" : "deathRate";
				else
					category = this.perturbIdx == 0 ? "edgeMean" : "edgeCV";
			}
			this.stats.increment(false, category);
		}

		// Specify the rate parameter
		double[] rates;
		if (this.level == Parameters.GENE_LEVEL)
			rates = pi.geneBDRates;
		else
			rates = (this.type == Parameters.BD_RATES) ? pi.domBDRates[this.domIndex] : pi.domEdgeRates[this.domIndex];

		rates[0] = this.ratesCache[0];
		rates[1] = this.ratesCache[1];

		this.ratesCache = null;
		this.perturbIdx = -1;
	}

	// Returns the "forward" probability density Q(x';x) for obtaining the new value
	// x' given the old value x.
	@Override
	public LogDouble getForwardDensity() {
		return this.forwardDensity;
	}

	// Returns the "backward" probability density Q(x;x') for obtaining the old
	// value x given the new value x'.
	@Override
	public LogDouble getBackwardDensity() {
		return this.backwardDensity;
	}

	// The ratio between the "backward" and "forward" proposal densities
	// respectively.
	// It will returns the ratio Q(x;x')/Q(x';x) for the old state x and the new
	// state x',
	@Override
	public LogDouble getDensityRatio() {
		return this.backwardDensity.divToNew(this.forwardDensity);
	}

	// This will return the index associated with perturbed rates
	public int getPerturbedRateIdx() {
		return this.perturbIdx;
	}

	// This method will return the name of acting proposer
	@Override
	public String getProposerName() {
		return this.name;
	}

	// This method will return the object of statistics associated with this
	// proposer
	@Override
	public ProposerStatistics getProposerStatistics() {
		return this.stats;
	}

	public void setProposalCV(double proposalCV) {
		this.proposalCV = proposalCV;
	}

	@Override
	public boolean isEnabled() {
		return this.isEnabled;
	}

	@Override
	public void setEnabled(boolean isActive) {
		this.isEnabled = isActive;
	}

	@Override
	public boolean hasValidProposal() {
		return (!this.backwardDensity.isZero() && !this.forwardDensity.isZero());
	}

}

/*
 * if (this.stats != null) { String category = (this.perturbIdx == 0) ?
 * "birthRate": "deathRate"; this.stats.increment(false, category); }
 */
// System.out.println("\n");
// if(this.level != Parameters.GENE_LEVEL && this.type != Parameters.BD_RATES)
// System.out.println("domain index =" + this.domIndex + " local cache = " +
// this.ratesCache[0] + "\t" + this.ratesCache[1] );
// if(this.level != Parameters.GENE_LEVEL && this.type != Parameters.BD_RATES)
// System.out.println("domain index =" + this.domIndex + " Global pi = " +
// pi.domEdgeRates[this.domIndex][0] + "\t" + pi.domEdgeRates[this.domIndex][1]
// );
