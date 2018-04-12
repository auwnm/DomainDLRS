package mcmc;

import java.util.ArrayList;

import proposers.Proposer;
import common.LogDouble;
import common.PRNG;

public class MetropolisHastingsAcceptor implements ProposalAcceptor {

	// Pseudo-random number generator
	private PRNG prng;

	// Constructor
	public MetropolisHastingsAcceptor(PRNG prng) {
		this.prng = prng;
	}

	// Returns true if a proposed state x' should be accepted according to the
	// Metropolis-Hastings sampling scheme
	@Override
	public boolean acceptProposedState(LogDouble proposedStateLikelihood, LogDouble oldStateLikelihood,
			ArrayList<Proposer> proposals) {
		LogDouble a = proposedStateLikelihood.divToNew(oldStateLikelihood);
		if (proposals != null) {
			for (Proposer prop : proposals) {
				if (!prop.hasValidProposal()) {
					return false;
				}
				a.mult(prop.getDensityRatio());
			}
		}
		return a.greaterThanOrEquals(new LogDouble(prng.nextDouble())); // Accounts also for case a >= 1.0.
	}

	@Override
	public boolean hasExhausted() {
		return false;
	}

}
