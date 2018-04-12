package mcmc;

import java.util.ArrayList;
import proposers.Proposer;
import common.LogDouble;

public interface ProposalAcceptor {

	// Returns true if a proposed state x' should be accepted, or false if the old
	// state x should be retained.
	public boolean acceptProposedState(LogDouble proposedStateLikelihood, LogDouble oldStateLikelihood,
			ArrayList<Proposer> proposals);

	// Return true only in Hill Climbing when it has exhausted maximum allowable
	// tries.
	public boolean hasExhausted();

}
