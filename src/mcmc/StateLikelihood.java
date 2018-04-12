package mcmc;

import common.LogDouble;

public class StateLikelihood {

	public LogDouble unnormalizedDensity;
	public int numberOfDomains;
	public LogDouble geneDLModel;
	public LogDouble domainDLRModel;
	public LogDouble[] domainSeqEvoModel;

	public StateLikelihood(int numberOfDomains) {
		this.geneDLModel = null;
		this.domainDLRModel = null;
		this.numberOfDomains = numberOfDomains;
		this.domainSeqEvoModel = new LogDouble[numberOfDomains];
	}

	public void updateOverall() {
		unnormalizedDensity = new LogDouble(1.0);
		unnormalizedDensity.mult(geneDLModel);
		unnormalizedDensity.mult(domainDLRModel);
		for (int r = 0; r < numberOfDomains; r++)
			unnormalizedDensity.mult(domainSeqEvoModel[r]);
	}

	public LogDouble getUnnormalizedDensity() {
		return unnormalizedDensity;
	}

	public String show() {

		StringBuilder sb = new StringBuilder();

		sb.append("[");
		sb.append(String.format(" %3.8f %3.8f", this.geneDLModel.getLogValue(), this.domainDLRModel.getLogValue()));
		for (int r = 0; r < this.numberOfDomains; r++)
			sb.append(String.format(" %3.8f", this.domainSeqEvoModel[r].getLogValue()));

		sb.append(String.format(" %3.8f", this.unnormalizedDensity.getLogValue()));
		sb.append("]");

		return sb.toString();

	}

	/*
	 * System.out.println("Gene DL Model density = " + this.geneDLModel);
	 * System.out.println("Domain DLR Model density = " + this.domainDLRModel);
	 * for(int r=0 ; r< numberOfDomains ; r++ ) System.out.println("Domain_" + (r+1)
	 * +" Seq Evo Model likelihood = " + this.domainSeqEvoModel[r]);
	 * System.out.println("Overall Likelihood = " + this.overall);
	 */

}
