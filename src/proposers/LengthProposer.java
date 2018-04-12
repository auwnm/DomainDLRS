package proposers;

import java.util.ArrayList;

import Main.Parameters;
import Main.TuningParameters;
import common.LogDouble;
import common.Mathematics;
import common.PRNG;
import common.RealInterval;

public class LengthProposer implements Proposer {

	public String name;
	public Parameters pi;
	public int domIndex;

	public PRNG prng;
	public double proposalCV;
	public RealInterval interval;
	public ProposerStatistics stats;
	public double[] lengthsSelectorWeights;

	public ArrayList<Integer> perturbIdx; // Indicies of the value going to be perturbed

	private LogDouble forwardDensity; // Forward probability density.
	private LogDouble backwardDensity; // Backward probability density.
	private boolean isEnabled; // On/off switch. By default it will be true.
	private boolean useRootBranch;

	private static double minBranchLength = TuningParameters.minBranchLength;

	public LengthProposer(Parameters pi, int domIndex, double proposalCV, PRNG prng, boolean useRootBranch) {
		this.name = "domBranchLengthProp_" + (domIndex + 1);
		this.pi = pi;
		this.domIndex = domIndex;
		this.prng = prng;
		this.proposalCV = proposalCV;
		this.isEnabled = true;
		this.interval = new RealInterval(minBranchLength, Double.POSITIVE_INFINITY, true, true);
		this.lengthsSelectorWeights = Mathematics.prepareCumulative(TuningParameters.lengthsSelectorWeights);
		this.stats = new ProposerStatistics(this.name);
		this.useRootBranch = useRootBranch;
	}

	// This method will perturb branch lengths w.r.t. its tree structure
	// nextInt = Returns a pseudorandom, uniformly distributed int value between 0
	// (inclusive) and the specified value (exclusive)
	@Override
	public boolean cacheAndPerturb() {

		// Specify the parameter
		double[] branchLength = this.pi.domainTree[this.domIndex].bl;

		// First cache the current state of length parameter(s)
		this.pi.domainTree[this.domIndex].cacheLengths();

		// Determine the number of lengths that will be perturbed
		double d = this.prng.nextDouble();
		int noOfLens = 1;
		while (d > this.lengthsSelectorWeights[noOfLens - 1]) {
			++noOfLens;
		}

		// Randomly pick branches for the perturbation
		this.perturbIdx = new ArrayList<Integer>();
		do {
			int idx = this.prng.nextInt(branchLength.length);

			if (this.pi.domainTree[this.domIndex].root.id == idx && !this.useRootBranch)
				continue;

			if (perturbIdx.contains(idx))
				continue;

			perturbIdx.add(idx);

		} while (perturbIdx.size() != noOfLens);

		// Perturbing the selected indicies
		LogDouble[] fbDensities = NormalProposer.perturbMultiValues(branchLength, perturbIdx, interval, this.proposalCV,
				prng, false);
		this.forwardDensity = fbDensities[0];
		this.backwardDensity = fbDensities[1];

		return true;
	}

	// Clears the cached (previous) state of perturbed parameters when a proposed
	// state has been accepted.
	@Override
	public void clearCache() {
		if (this.stats != null) {
			String category = "" + this.perturbIdx.size() + " perturbed length sub-parameters";
			this.stats.increment(true, category);
		}
		this.pi.domainTree[this.domIndex].clearLengthsCache();
		this.perturbIdx = null;
	}

	// Restores the cached (previous) state of perturbed parameters when a proposed
	// state has been rejected.
	@Override
	public void restoreCache() {
		if (this.stats != null) {
			String category = "" + this.perturbIdx.size() + " perturbed length sub-parameters";
			this.stats.increment(false, category);
		}

		this.pi.domainTree[this.domIndex].restoreLengths();
		this.perturbIdx = null;
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

	// This will return the index associated with perturbed rates
	public ArrayList<Integer> getPerturbedIdx() {
		return this.perturbIdx;
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

// Pair<double [],LogDouble[]> perturbedValues =
// NormalProposer.perturbMultiValues(branchLength, perturbIdx, interval,
// this.proposalCV, prng);
// if(perturbedValues == null)
// return false;
// Changing the state of parameter
// for(int i=0 ; i<perturbIdx.size() ; i++)
// this.lengths[ perturbIdx.get(i) ] = perturbedValues.first[perturbIdx.get(i)];
// this.forwardDensity = perturbedValues.second[0];
// this.backwardDensity = perturbedValues.second[1];
/*
 * public void setLengths(double [] rates) { this.lengths = rates; }
 */