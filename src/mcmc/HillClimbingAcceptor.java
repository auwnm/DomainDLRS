package mcmc;

import java.util.ArrayList;
import proposers.Proposer;
import common.LogDouble;

public class HillClimbingAcceptor implements ProposalAcceptor {

	private long maxTries; // Max number of non-improvements.
	private long noOfTries; // Number of straight non-improvements.
	private boolean exhausted = false; // Tried all allowable maximum tries
	private double precision = 1.0 + 1e-32; // Precision for improvement comparison.

	public HillClimbingAcceptor(long maxTries) {
		this.noOfTries = 0;
		this.maxTries = maxTries;
		this.exhausted = false;
	}

	// Returns true if a proposed state x' has a higher likelihood than the old
	// state x.
	// Does not require any proposal information nor a pseudo-random number
	// generator.
	@Override
	public boolean acceptProposedState(LogDouble proposedStateLikelihood, LogDouble oldStateLikelihood,
			ArrayList<Proposer> proposals) {

		// If we do have proposal information, just verify it's OK.
		if (proposals != null) {
			for (Proposer prop : proposals) {
				if (!prop.hasValidProposal()) {
					return false;
				}
			}
		}

		boolean isImproved = proposedStateLikelihood.greaterThan(oldStateLikelihood.multToNew(this.precision));
		if (isImproved) {
			this.noOfTries = 0;
		} else {
			this.noOfTries++;
		}

		if (this.noOfTries > this.maxTries)
			this.exhausted = true;

		return isImproved;
	}

	@Override
	public boolean hasExhausted() {
		return this.exhausted;
	}

}
