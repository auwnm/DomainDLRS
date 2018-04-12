package Main;

import phylogeny.Mapping;
import phylogeny.Tree;
import seqEvolution.MSAData;
import seqEvolution.SubstitutionMatrix;

public class Parameters {

	public static int GENE_LEVEL = 0;
	public static int DOMAIN_LEVEL = 1;

	public static int BD_RATES = 0;
	public static int EDGE_RATES = 1;
	public static int LENGTHS = 2;
	public static int TOPOLOGY = 3;

	public int numberOfDomains;
	public SubstitutionMatrix Q;
	public Tree speciesTree;
	public Tree geneTree;
	public Tree[] domainTree;
	public MSAData[] domMSA;

	public double[] geneBDRates;
	public double[][] domBDRates;
	public double[][] domEdgeRates;

	public Mapping gene2species;
	public Mapping[] domain2genes;

	public Parameters(int numberOfDomains) {
		this.speciesTree = null;
		this.geneTree = null;
		this.gene2species = null;
		this.domMSA = null;
		this.Q = null;
		this.numberOfDomains = numberOfDomains;
		this.geneBDRates = new double[2];
		this.domBDRates = new double[numberOfDomains][2];
		this.domEdgeRates = new double[numberOfDomains][2];
		this.domainTree = new Tree[numberOfDomains];
		this.domain2genes = new Mapping[numberOfDomains];
	}

	// This method will show the state of parameters.
	public void showState(int level, int paramType, int domIndex) {

		// Construct stateStr to append
		StringBuilder stateStr = new StringBuilder();

		if (level == GENE_LEVEL) {
			if (paramType == BD_RATES)
				stateStr.append(String.format("[ %f %f ]", this.geneBDRates[0], this.geneBDRates[1]));
			if (paramType == TOPOLOGY) {
				stateStr.append("parray = [");
				for (int i = 0; i < this.geneTree.nnodes; i++)
					stateStr.append(this.geneTree.parray[i] + ",");
				stateStr.append("]");
			}

		} else {
			if (paramType == BD_RATES) {
				stateStr.append(
						String.format("\t [%f %f]", this.domBDRates[domIndex][0], this.domBDRates[domIndex][1]));
			} else if (paramType == EDGE_RATES) {
				stateStr.append(
						String.format("\t [%f %f]", this.domEdgeRates[domIndex][0], this.domEdgeRates[domIndex][1]));
			} else if (paramType == TOPOLOGY) {
				stateStr.append("[");
				for (int i = 0; i < this.domainTree[domIndex].nnodes; i++)
					stateStr.append(this.domainTree[domIndex].parray[i] + ",");
				stateStr.append("]\t");

			} else if (paramType == LENGTHS) {
				stateStr.append("[");
				for (int i = 0; i < this.domainTree[domIndex].nnodes; i++)
					stateStr.append(this.domainTree[domIndex].bl[i] + ",");
				stateStr.append("]\t");
			}

		}
		System.out.println(stateStr);
	}

}

/*
 * // Updating Parameter ////////////////////////////////// public void
 * updateGeneBDRates(double [] rates) { for(int i=0 ; i< rates.length ; i++)
 * this.geneBDRates[i] = rates[i]; }
 * 
 * public void updateDomainBDRates(int domIdx,double [] rates) { for(int i=0 ;
 * i< rates.length ; i++) this.domBDRates[domIdx][i] = rates[i]; }
 * 
 * public void updateDomainEdgeRates(int domIdx,double [] rates) { for(int i=0 ;
 * i< rates.length ; i++) this.domEdgeRates[domIdx][i] = rates[i]; }
 * 
 * 
 * // Parameter-wise cache ////////////////////////////////// private Tree []
 * treeCache; private double [][] rateCache; private ArrayList<double []>
 * lengthCache; this.treeCache = null; this.lengthCache = null; this.rateCache =
 * null;
 * 
 * // Tree Parameter public boolean cacheTrees() { this.treeCache = new
 * Tree[this.numberOfDomains + 1]; int r=0; for(;r<this.numberOfDomains;r++)
 * this.treeCache[r] = this.domainTree[r].copy(); this.treeCache[r] =
 * this.geneTree.copy(); return true; } public boolean restoreTrees() { int r=0;
 * for(;r<this.numberOfDomains;r++) this.domainTree[r] = this.treeCache[r];
 * this.geneTree = this.treeCache[r]; return true; } public boolean
 * clearTreeCache() { this.treeCache = null; return true; }
 * 
 * 
 * 
 * // Rates Parameter public boolean cacheRates() { int noOfRateParameters = 4;
 * this.rateCache = new double [this.numberOfDomains + 1][noOfRateParameters];
 * int r=0; for(;r<this.numberOfDomains;r++) { this.rateCache[r][0] =
 * this.domBDRates[r][0]; this.rateCache[r][1] = this.domBDRates[r][1];
 * this.rateCache[r][2] = this.domEdgeRates[r][0]; this.rateCache[r][3] =
 * this.domEdgeRates[r][1]; } this.rateCache[r][0] = this.geneBDRates[0];
 * this.rateCache[r][1] = this.geneBDRates[1];
 * 
 * return true; } public boolean restoreRates() { int r=0;
 * for(;r<this.numberOfDomains;r++) { this.domBDRates[r][0] =
 * this.rateCache[r][0]; this.domBDRates[r][1] = this.rateCache[r][1];
 * this.domEdgeRates[r][0] = this.rateCache[r][2]; this.domEdgeRates[r][1] =
 * this.rateCache[r][3]; } this.geneBDRates[0] = this.rateCache[r][0];
 * this.geneBDRates[1] = this.rateCache[r][1]; return true; } public boolean
 * clearRateCache() { this.rateCache = null; return true; }
 * 
 * 
 * 
 * // Length Parameter public boolean cacheLengths() { lengthCache = new
 * ArrayList<double[]>(); for(int r=0 ; r<this.numberOfDomains ; r++)
 * lengthCache.add(Arrays.copyOf(this.domainTree[r].bl,
 * this.domainTree[r].bl.length));
 * 
 * return true; } public boolean restoreLengths() { for(int r=0 ;
 * r<this.numberOfDomains ; r++) this.domainTree[r].bl = lengthCache.get(r);
 * return true; } public boolean clearLengthCache() { this.lengthCache = null;
 * return true; }
 * 
 * 
 * // Complete Cache ////////////////////////////////// Parameters cache;
 * this.cache = null; public boolean cache() {
 * 
 * cache = new Parameters(this.numberOfDomains);
 * 
 * cache.numberOfDomains = this.numberOfDomains; cache.Q = this.Q; cache.domMSA
 * = this.domMSA; cache.gene2species = this.gene2species; cache.domain2genes =
 * this.domain2genes; cache.speciesTree = this.speciesTree;
 * 
 * 
 * cache.geneTree = geneTree.copy(); for(int r=0 ;r<this.numberOfDomains ; r++)
 * cache.domainTree[r] = this.domainTree[r].copy();
 * 
 * cache.geneBDRates = Arrays.copyOf(this.geneBDRates, this.geneBDRates.length);
 * for(int r=0 ;r<this.numberOfDomains ; r++) { cache.domBDRates[r] =
 * Arrays.copyOf(this.domBDRates[r], 2); cache.domEdgeRates[r] =
 * Arrays.copyOf(this.domEdgeRates[r], 2); }
 * 
 * cache.cache = null;
 * 
 * return true; }
 * 
 * public void restore() {
 * 
 * this.numberOfDomains = cache.numberOfDomains; this.Q = cache.Q; this.domMSA =
 * cache.domMSA; this.gene2species = cache.gene2species; this.domain2genes =
 * cache.domain2genes; this.speciesTree = cache.speciesTree;
 * 
 * 
 * this.geneTree = cache.geneTree; for(int r=0 ;r<cache.numberOfDomains ; r++)
 * this.domainTree[r] = cache.domainTree[r];
 * 
 * this.geneBDRates = Arrays.copyOf(cache.geneBDRates,
 * cache.geneBDRates.length); for(int r=0 ;r<this.numberOfDomains ; r++) {
 * this.domBDRates[r] = Arrays.copyOf(cache.domBDRates[r], 2);
 * this.domEdgeRates[r] = Arrays.copyOf(cache.domEdgeRates[r], 2); }
 * 
 * this.cache = null;
 * 
 * }
 * 
 * public void clearCache() { this.cache = null; }
 */
