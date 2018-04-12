package proposers;

import java.util.ArrayList;
import java.util.HashMap;

import Main.Parameters;
import Main.TuningParameters;

import common.Mathematics;
import common.PRNG;

public class ProposerSelector {

	public int numberOfDomains;
	private double domainLevelWeight;
	private double[] noProposersWeights;
	private double[] domParametersWeights;
	private double[] geneParametersWeights;
	private double[] domainSelectorWeights;

	public int level;

	private PRNG prng;
	public ArrayList<Proposer> geneLevelProposers;
	public HashMap<Integer, ArrayList<Proposer>> domLevelProposers;

	public ProposerSelector(int numberOfDomains, PRNG prng) {
		this.prng = prng;
		this.numberOfDomains = numberOfDomains;
		this.domainLevelWeight = TuningParameters.domainLevelWeight;
		this.domainSelectorWeights = Mathematics.prepareCumulative(TuningParameters.domainSelectorWeights);
		this.noProposersWeights = Mathematics.prepareCumulative(TuningParameters.noProposersWeights);
		this.domParametersWeights = Mathematics.prepareCumulative(TuningParameters.domParametersWeights);
		this.geneParametersWeights = Mathematics.prepareCumulative(TuningParameters.geneParametersWeights);
		this.geneLevelProposers = new ArrayList<Proposer>();
		this.domLevelProposers = new HashMap<Integer, ArrayList<Proposer>>();
	}

	// Select the level of perturbation i.e. gene level or at domain level
	private int SelectLevel() {
		if (this.prng.nextDouble() < this.domainLevelWeight)
			return Parameters.DOMAIN_LEVEL;
		else
			return Parameters.GENE_LEVEL;
	}

	// Uniform Selection of domain tree
	private int SelectDomain() {

		if (this.numberOfDomains == 1)
			return 0;

		double rv = this.prng.nextDouble();
		int domIdx = 0;
		while (rv > this.domainSelectorWeights[domIdx]) {
			++domIdx;
		}
		return domIdx;

		// return this.prng.nextInt(this.numberOfDomains);
	}

	// This method will propose the list of perturbers
	public void SelectProposers(ArrayList<Proposer> proposalList) {

		// First determine the level of perturbation i.e. at gene level or at domain
		// level
		level = SelectLevel();

		if (level == Parameters.GENE_LEVEL) {
			int noOfProps = 1;
			do {
				double rv = this.prng.nextDouble();
				int paramIdx = 0;
				while (rv > this.geneParametersWeights[paramIdx]) {
					++paramIdx;
				}
				Proposer p = this.geneLevelProposers.get(paramIdx);
				if (!proposalList.contains(p) && p.isEnabled())
					proposalList.add(p);
			} while (proposalList.size() != noOfProps);

		} else {

			// Pick domain tree uniformly
			// Select only one domain (dimension) at each perturbation
			int domIdx = SelectDomain();

			// Check the number of active proposers for this domain
			int noOfActiveProposers = 0;
			for (int a = 0; a < this.domLevelProposers.get(domIdx).size(); a++) {
				Proposer p = this.domLevelProposers.get(domIdx).get(a);
				if (p.isEnabled())
					noOfActiveProposers++;
			}

			// Determine the number of proposers that will take part in perturbation at
			// domain level
			double d = this.prng.nextDouble();
			int noOfProps = 1;
			while (d > this.noProposersWeights[noOfProps - 1]) {
				++noOfProps;
			}

			// In case, if it is greater than the number of active proposers
			if (noOfProps > noOfActiveProposers)
				noOfProps = noOfActiveProposers;

			// Loop until desired number of proposers collected ...
			do {
				double rv = this.prng.nextDouble();
				int paramIdx = 0;
				while (rv > this.domParametersWeights[paramIdx]) {
					++paramIdx;
				}
				Proposer p = this.domLevelProposers.get(domIdx).get(paramIdx);

				if (!proposalList.contains(p) && p.isEnabled()) {
					proposalList.add(p);
				}

			} while (proposalList.size() != noOfProps);

		}

		return;
	}

}

/*
 * // This method will propose the list of perturbers public void
 * SelectProposers(ArrayList<Proposer> proposalList) {
 * 
 * // Determine the number of propopsers that will take part in perturbation
 * double d = this.prng.nextDouble(); int noOfProps = 1; while (d >
 * this.noProposersWeights[noOfProps-1]) { ++noOfProps; } noOfProps *=
 * numberOfDomains;
 * 
 * // Loop until desired number of proposers collected ... do {
 * 
 * // First determine the level of perturbation i.e. at gene level or at domain
 * level level = SelectLevel(); if(level == Parameters.GENE_LEVEL) { double rv =
 * this.prng.nextDouble(); int paramIdx = 0; while ( rv >
 * this.geneParametersWeights[paramIdx] ) { ++paramIdx; } Proposer p =
 * this.geneLevelProposers.get(paramIdx); if( !proposalList.contains(p) &&
 * p.isEnabled() ) proposalList.add(p); } else {
 * 
 * // Pick domain tree uniformly int domIdx = SelectDomain(); double rv =
 * this.prng.nextDouble(); int paramIdx = 0; while ( rv >
 * this.domParametersWeights[paramIdx] ) { ++paramIdx; } Proposer p =
 * this.domLevelProposers.get(domIdx).get(paramIdx);
 * if(!proposalList.contains(p) && p.isEnabled() ) proposalList.add(p); }
 * 
 * } while (proposalList.size() != noOfProps);
 * 
 * return ; }
 */
