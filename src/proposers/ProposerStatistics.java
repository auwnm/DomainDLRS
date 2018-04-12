package proposers;

import java.io.IOException;
import java.util.Map.Entry;
import java.util.TreeMap;

import common.Log;

/**
 * Base class for statistics of <code>Proposer</code> objects (and potentially
 * others) with respect to acceptance and rejection. If desired, one may
 * increment the acceptance/rejections counter by also providing a string key
 * for a specific sub-category that should be monitored.
 * <p/>
 * See also <code>FineProposerStatistics</code>.
 * 
 * @author Joel Sjöstrand.
 */

public class ProposerStatistics {

	private String name; // Associated proposer name
	private int noOfAccepted; // Overall number of accepted proposals
	private int noOfRejected; // Overall number of rejected proposals.
	private TreeMap<String, int[]> accRejByKey; // Hash for key-specific acceptance/rejections.

	// Constructor starting with associated proposer name
	public ProposerStatistics(String name) {
		this.name = name;
		this.noOfAccepted = 0;
		this.noOfRejected = 0;
		this.accRejByKey = new TreeMap<String, int[]>();
	}

	// Returns the number of times the associated Proposer has suggested new states.
	public int getNoOfProposals() {
		return (this.noOfAccepted + this.noOfRejected);
	}

	// Returns the number of times the associated <code>Proposer</code> has accepted
	// new states.
	public int getNoOfAcceptedProposals() {
		return this.noOfAccepted;
	}

	// Returns the number of times the associated <code>Proposer</code> has rejected
	// new states.
	public int getNoOfRejectedProposals() {
		return this.noOfRejected;
	}

	// Returns the number of accepted proposals divided by the total number of
	// proposals.
	public double getAcceptanceRatio() {
		return (this.noOfAccepted / (double) (this.noOfAccepted + this.noOfRejected));
	}

	// Returns the number of rejected proposals divided by the total number of
	// proposals.
	public double getRejectionRatio() {
		return (this.noOfRejected / (double) (this.noOfAccepted + this.noOfRejected));
	}

	// Adds a proposal outcome.
	// If new state was accepted; wasAccepted = true otherwise false if rejected.
	public void increment(boolean wasAccepted) {
		if (wasAccepted) {
			++this.noOfAccepted;
		} else {
			++this.noOfRejected;
		}
	}

	// Adds a proposal outcome.
	// If new state was accepted; wasAccepted = true otherwise false if rejected.
	// Category means specific sub-category to which the acceptance/rejection
	// belongs.
	public void increment(boolean wasAccepted, String category) {
		this.increment(wasAccepted);
		int[] cat = this.accRejByKey.get(category);
		if (cat == null) {
			cat = new int[2];
			cat[wasAccepted ? 0 : 1] = 1;
			this.accRejByKey.put(category, cat);
		} else {
			cat[wasAccepted ? 0 : 1] += 1;
		}
	}

	// Log the statistics of this proposer
	public boolean writeProposerStats(Log log) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("Proposer stats:\t" + this.name + "\n");
		sb.append("Acceptance ratio: ").append(this.noOfAccepted).append(" / ").append(this.getNoOfProposals())
				.append(" = ").append(this.getAcceptanceRatio()).append("\n");
		if (!this.accRejByKey.isEmpty()) {
			sb.append("Acceptance ratios for sub-categories:\n");
			for (Entry<String, int[]> kv : this.accRejByKey.entrySet()) {
				int acc = kv.getValue()[0];
				int rej = kv.getValue()[1];
				sb.append('\t').append(kv.getKey()).append('\t').append(acc).append(" / ").append(acc + rej)
						.append(" = ").append(acc / (double) (acc + rej)).append("\n");
			}
		}
		if (log != null)
			log.write(sb.toString());
		return true;
	}

}
