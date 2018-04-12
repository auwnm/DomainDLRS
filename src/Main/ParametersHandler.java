package Main;

import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import mcmc.RealParameterUniformPrior;

import common.Log;
import common.LogDouble;
import common.PRNG;
import common.Pair;

import dlModel.BirthDeathProbs;
import dlModel.DLModel;

import phylogeny.Mapping;
import phylogeny.NewickReader;
import phylogeny.Node;
import phylogeny.Reconciliation;

import phylogeny.Tree;

import rateModel.RateModel;
import seqEvolution.FastaRW;
import seqEvolution.MSAData;
import seqEvolution.SMD;
import seqEvolution.SequenceType;
import seqEvolution.Sequences;
import seqEvolution.SubstitutionMatrix;

import treeConstructor.DistanceMatrix;
import treeConstructor.NJConstructor;

public class ParametersHandler {

	public static Log errorLog = null;

	// Set Error log for parameter handler
	public static void setErrorLog(Log errorLog) {
		ParametersHandler.errorLog = errorLog;
	}

	// Initializing Substitution model
	public static boolean initSubstitutionModel(Parameters pi, String substitutionModel) throws IOException {

		SubstitutionMatrix Q = SMD.create(substitutionModel);
		if (Q == null && errorLog != null) {
			errorLog.write("\nParameterHandler Error: Specified substitution model is not included.");
			return false;
		}

		pi.Q = Q;
		return true;
	}

	// Reading gene mapping file
	public static boolean readGeneMappingFile(Parameters pi, String path, String geneMapFile) throws IOException {

		Mapping gene2species = new Mapping();
		boolean success = gene2species.readMapping(path + geneMapFile);
		if (!success || gene2species.map == null) {
			errorLog.write("\nParameterHandler Error: Cannot read mapping file " + geneMapFile);
			return false;
		}

		pi.gene2species = gene2species;
		return true;
	}

	// Reading domain mapping files
	public static boolean readDomainMappingFiles(Parameters pi, String path, String[] mappingFiles) throws IOException {

		int noOfMappingFiles = mappingFiles.length;
		Mapping[] mappings = new Mapping[noOfMappingFiles];
		for (int i = 0; i < noOfMappingFiles; i++) {
			mappings[i] = new Mapping();
			boolean success = mappings[i].readMapping(path + mappingFiles[i]);
			if (!success || mappings[i].map == null) {
				errorLog.write("\nParameterHandler Error: Cannot read mapping file " + mappingFiles[i]);
				return false;
			}
		}

		pi.domain2genes = mappings;
		return true;
	}

	// Reading alignment files
	public static boolean readAlignmentFiles(Parameters pi, String path, String[] alignmentFiles) throws IOException {

		// Reading domain alignment files one by one
		int noOfAlignmentFiles = alignmentFiles.length;
		Sequences[] alignments = new Sequences[noOfAlignmentFiles];
		for (int i = 0; i < noOfAlignmentFiles; i++) {
			alignments[i] = FastaRW.readFastaAlign(path + alignmentFiles[i]);
			if (alignments[i] == null) {
				errorLog.write("\nParameterHandler Error: Cannot read sequence file " + alignmentFiles[i]);
				return false;
			}
		}
		// Making MSA objects from these alignments
		SubstitutionMatrix Q = pi.Q;
		MSAData[] domainMSA = new MSAData[noOfAlignmentFiles];
		for (int i = 0; i < noOfAlignmentFiles; i++)
			domainMSA[i] = new MSAData(Q.getSequenceType(), alignments[i]);

		pi.domMSA = domainMSA;
		return true;
	}

	// Initialize the gene tree by given tree or by using NJ on random domain based
	// MSA or by random topology
	public static boolean initGeneTree(Parameters pi, String geneTree, PRNG prng, boolean neighborJoining)
			throws IOException {
		if (geneTree != null) {
			pi.geneTree = NewickReader.readNewickTreeStr(geneTree);
			pi.geneTree.setName("geneTree");
			pi.geneTree.setDepths(pi.geneTree.root, 0);
			return true;
		}
		if (neighborJoining) {
			// Constructing gene tree by NJ algorithm on random MSA
			MSAData randomGeneMSA = getGeneMSA(pi.domMSA, pi.domain2genes, pi.gene2species, pi.Q.getSequenceType(),
					prng);
			pi.geneTree = constructTree(randomGeneMSA, pi.Q, prng);
			pi.geneTree.setName("geneTree");
			pi.geneTree.setDepths(pi.geneTree.root, 0);

		} else {
			// Random topology construction of gene tree
			String[] leavesName = new String[pi.gene2species.map.size()];
			int idx = 0;
			for (String key : pi.gene2species.map.keySet())
				leavesName[idx++] = key;
			pi.geneTree = NJConstructor.RandomTree(leavesName, prng, 0);
			pi.geneTree.setName("geneTree");
			pi.geneTree.setDepths(pi.geneTree.root, 0);
		}
		return true;
	}

	// Initialize the domain tree
	public static boolean initDomainTrees(Parameters pi, PRNG prng, boolean neighborJoining) throws IOException {
		if (neighborJoining) {
			// Constructing domain trees with their lengths by using NJ algorithm
			for (int i = 0; i < pi.numberOfDomains; i++) {
				pi.domainTree[i] = ParametersHandler.constructTree(pi.domMSA[i], pi.Q, prng);
				pi.domainTree[i].setName("domainTree" + (i + 1));
			}
		} else {
			// Random topology construction of domain trees
			for (int i = 0; i < pi.numberOfDomains; i++) {
				pi.domainTree[i] = ParametersHandler.constructRandomTree(pi.domMSA[i], prng);
				pi.domainTree[i].setName("domainTree" + (i + 1));
			}
		}
		return true;
	}

	// Read the time species tree
	public static boolean readSpeciesTree(Parameters pi, String path, String speciesTreeFile) throws IOException {
		Tree speciesTree = NewickReader.readNewickTreeFile(path + speciesTreeFile);
		speciesTree.setName("speciesTree");
		speciesTree.setDepths(speciesTree.root, 0);
		speciesTree.setTimesFromOriginalLengths();
		speciesTree.setNormalizedTimes();
		pi.speciesTree = speciesTree;
		return true;
	}

	// Returns a PRNG. If no seed is found, a random seed is used.
	public static PRNG getPRNG(String seed) {
		return (seed == null ? new PRNG() : new PRNG(new BigInteger(seed)));
	}

	// Construct random tree
	public static Tree constructRandomTree(MSAData geneMSA, PRNG prng) {
		double lengthScale = 1.0;
		Tree tree = NJConstructor.RandomTree(geneMSA.getSequencesIdentifiers(), prng, lengthScale);
		return tree;
	}

	// Construct tree based on Neighbour Joining Technique
	public static Tree constructTree(MSAData geneMSA, SubstitutionMatrix Q, PRNG prng) {

		double[][] distances = DistanceMatrix.compute(geneMSA, Q, DistanceMatrix.JC69);
		Tree tree = NJConstructor.RasmussenNJ(distances, geneMSA.getSequencesIdentifiers());

		// To handle the negative and zero branch lengths
		double[] lengths = tree.bl;
		for (int i = 0; i < tree.nnodes; i++) {
			if (lengths[i] <= TuningParameters.minBranchLength)
				lengths[i] = prng.nextDouble();
		}

		return tree;
	}

	// This method will initialize the paramters with actual simulation parameters.
	public static void InitializeParametersWithActualEstimates(Parameters pi) {

		// Scaling factor for adjusting birth death rates according to time
		// normalization
		double scalingFactor = 1.0; // 142.7;

		// Birth death parameters
		pi.geneBDRates[0] = 0.012 * scalingFactor;
		pi.geneBDRates[1] = 0.010 * scalingFactor;
		pi.domBDRates[0][0] = 0.012125323499552323 * scalingFactor;
		pi.domBDRates[0][1] = 0.010687097698647693 * scalingFactor;
		pi.domBDRates[1][0] = 0.012342038623025768 * scalingFactor;
		pi.domBDRates[1][1] = 0.010796581085432334 * scalingFactor;

		// Edge rate parameters
		double k = 0.49;
		double theta = 1.1;
		pi.domEdgeRates[0][0] = k * theta;
		pi.domEdgeRates[0][1] = Math.sqrt(k) / k;
		pi.domEdgeRates[1][0] = k * theta;
		pi.domEdgeRates[1][1] = Math.sqrt(k) / k;

	}

	// MPR based rate parameters initialization
	public static void initRates(Parameters pi, RealParameterUniformPrior edgeRatePriors, PRNG prng,
			boolean useRootEdge, Log log) throws IOException {

		// Compute MPR reconciliations between domain,gene and species tree
		Reconciliation geneRecon = new Reconciliation(pi.geneTree, pi.speciesTree, pi.gene2species);
		Reconciliation[] domainRecon = new Reconciliation[pi.numberOfDomains];
		for (int i = 0; i < pi.numberOfDomains; i++)
			domainRecon[i] = new Reconciliation(pi.domainTree[i], pi.geneTree, pi.domain2genes[i]);

		// Estimate rates based on these mpr reconciliations
		ParametersHandler.estimateBDRates(pi, geneRecon, domainRecon, log);
		ParametersHandler.estimateEdgeRates(pi, geneRecon, domainRecon, prng, useRootEdge, log);

		// If initial rate estimates are violating the priors (Only for edge rate prior)
		LogDouble meanPrior, cvPrior;
		double a = edgeRatePriors.interval.getLowerBound();
		double b = edgeRatePriors.interval.getUpperBound();
		for (int r = 0; r < pi.numberOfDomains; r++) {
			meanPrior = edgeRatePriors.getPriorProbability(pi.domEdgeRates[r][0]);
			cvPrior = edgeRatePriors.getPriorProbability(pi.domEdgeRates[r][1]);
			if (meanPrior.isZero() || cvPrior.isZero()) {
				pi.domEdgeRates[r][0] = a + prng.nextDouble() * b;
				pi.domEdgeRates[r][1] = a + prng.nextDouble() * b;
			}
		}

		// Log the intialized rates
		log.newLine();
		StringBuilder str = new StringBuilder();
		str.append("Initial rate Estimates:\n");

		str.append(String.format("%-50s", "Estimated gene birth death rates"));
		str.append("  : " + pi.geneBDRates[0] + " \t " + pi.geneBDRates[1]);
		str.append("\n");

		String BDStr = String.format("%-50s", "Estimated doamin BD rates for domain tree ");
		for (int i = 0; i < pi.numberOfDomains; i++)
			str.append(BDStr + (i + 1) + " : " + pi.domBDRates[i][0] + " \t " + pi.domBDRates[i][1] + "\n");

		String EdgeStr = String.format("%-50s", "Estimated doamin edge rates for domain tree ");
		for (int i = 0; i < pi.numberOfDomains; i++)
			str.append(EdgeStr + (i + 1) + " : " + pi.domEdgeRates[i][0] + " \t " + pi.domEdgeRates[i][1] + "\n");
		str.append("\n");

		log.write(str.toString());
		log.newLine();
	}

	// MPR based birth death rate parameter estimation
	public static void estimateBDRates(Parameters pi, Reconciliation geneRecon, Reconciliation[] domainRecon, Log log)
			throws IOException {

		StringBuilder logStr = new StringBuilder();
		int noLosses;
		int noDuplications;

		// Total arc time of host (species) tree (including stem arc)
		double totime = pi.speciesTree.getTotalTime();
		if (Double.isNaN(totime))
			logStr.append("\nInitialisation Error: Total time is NAN in birth death rate parameters estimation.");

		// Gene birth death rate parameter estimation
		noDuplications = 0;
		for (Node gnode : geneRecon.guest.nodes) {
			if (geneRecon.eventMap[gnode.id] == Reconciliation.EVENT_DUPL)
				noDuplications++;
		}

		noLosses = geneRecon.noImpliedNodes;
		if (noDuplications == 0) {
			pi.geneBDRates[0] = (pi.geneTree.getNumberOfLeaves() - 1) / totime;
			if (noLosses == 0)
				pi.geneBDRates[1] = pi.geneBDRates[0];
			else
				pi.geneBDRates[1] = noLosses / totime;
		} else {
			pi.geneBDRates[0] = noDuplications / totime;
			if (noLosses == 0)
				pi.geneBDRates[1] = pi.geneBDRates[0];
			else
				pi.geneBDRates[1] = noLosses / totime;

		}
		logStr.append("\nNumber of gene duplications nodes = " + noDuplications);
		logStr.append("\nNumber of gene implied loss nodes = " + noLosses);

		// Domain birth death rate parameter estimation
		for (int r = 0; r < pi.numberOfDomains; r++) {

			noDuplications = 0;
			for (Node dnode : domainRecon[r].guest.nodes)
				if (domainRecon[r].eventMap[dnode.id] == Reconciliation.EVENT_DUPL)
					noDuplications++;

			noLosses = domainRecon[r].noImpliedNodes;
			if (noDuplications == 0) {
				pi.domBDRates[r][0] = (pi.domainTree[r].getNumberOfLeaves() - 1) / totime;
				if (noLosses == 0)
					pi.domBDRates[r][1] = pi.domBDRates[r][0];
				else
					pi.domBDRates[r][1] = noLosses / totime;
			} else {
				pi.domBDRates[r][0] = noDuplications / totime;
				if (noLosses == 0)
					pi.domBDRates[r][1] = pi.domBDRates[r][0];
				else
					pi.domBDRates[r][1] = noLosses / totime;

			}
			logStr.append("\nNumber of duplications nodes in " + pi.domainTree[r].name + " = " + noDuplications);
			logStr.append("\nNumber of implied nodes in" + pi.domainTree[r].name + " = " + noLosses);
		}
		if (log != null && logStr.length() != 0)
			log.write(logStr.toString() + "\n");
	}

	// Estimation of rate parameters based on lengths and time sampled from birth
	// death process
	public static void estimateEdgeRates(Parameters pi, Reconciliation geneRecon, Reconciliation[] domainRecon,
			PRNG prng, boolean useRootEdge, Log log) throws IOException {

		int totalIterations = 1000;
		StringBuilder logStr = new StringBuilder();

		// Getting into likelihood business
		LogDouble geneDLModelLikelihood, domainDLModelLikelihood, rateModelLikelihoods;

		// Initialize rate parameters with random values
		for (int r = 0; r < pi.numberOfDomains; r++) {
			pi.domEdgeRates[r][0] = Math.random();
			pi.domEdgeRates[r][1] = Math.random();
		}

		// Initialize the DL model
		double[] speciesExtinctionProb = BirthDeathProbs.calcExtinctionProb(pi.speciesTree, pi.geneBDRates);
		DLModel geneDLModel = new DLModel(geneRecon, pi.geneBDRates, speciesExtinctionProb, useRootEdge, prng, log);
		geneDLModelLikelihood = geneDLModel.birthDeathTreePrior();
		double[] extinctionProb;
		DLModel[] domDLModel = new DLModel[pi.numberOfDomains];
		for (int r = 0; r < pi.numberOfDomains; r++)
			domDLModel[r] = new DLModel(domainRecon[r], pi.domBDRates[r], null, useRootEdge, prng, log);

		// Initialize the rate model
		RateModel[] rateModel = new RateModel[pi.numberOfDomains];
		for (int i = 0; i < pi.numberOfDomains; i++)
			rateModel[i] = new RateModel(pi.domEdgeRates[i], useRootEdge);

		double[] r8;
		double length, time;
		LogDouble newLikelihood, maxLikelihood;
		maxLikelihood = new LogDouble(0.0);
		for (int n = 0; n < totalIterations; n++) {

			newLikelihood = new LogDouble(1.0);
			newLikelihood.mult(geneDLModelLikelihood);
			geneDLModel.sampleRealisation();
			for (int r = 0; r < pi.numberOfDomains; r++) {
				extinctionProb = BirthDeathProbs.calcExtinctionProb(pi.geneTree, pi.domBDRates[r]);
				domDLModel[r].updateDLModel(extinctionProb);
				domainDLModelLikelihood = domDLModel[r].birthDeathTreePrior();
				domDLModel[r].sampleRealisation();
				rateModel[r].update(domainRecon[r].guest);
				rateModelLikelihoods = rateModel[r].getLogLikelihood();

				newLikelihood.mult(domainDLModelLikelihood);
				newLikelihood.mult(rateModelLikelihoods);
			}

			// Pick time sample which has maximum likelihood
			if (newLikelihood.greaterThan(maxLikelihood)) {
				StringBuilder ratesStr = new StringBuilder();
				for (int r = 0; r < pi.numberOfDomains; r++) {
					r8 = new double[pi.domainTree[r].nnodes];
					for (Node node : pi.domainTree[r].nodes) {
						length = pi.domainTree[r].bl[node.id];
						if (node.isRoot())
							time = pi.domainTree[r].getPeakTime() - pi.domainTree[r].vt[node.id];
						else
							time = pi.domainTree[r].vt[node.parent.id] - pi.domainTree[r].vt[node.id];
						r8[node.id] = length / time;
					}

					// Update the rates based on this new sample
					Pair<Double, Double> rates = getMLERates(r8);
					pi.domEdgeRates[r][0] = rates.first;
					pi.domEdgeRates[r][1] = rates.second;

					ratesStr.append(String.format("\t%1.5f", rates.first));
					ratesStr.append(String.format("\t%1.5f", rates.second));
				}
				logStr.append("\ni = " + String.format("%-6d", n) + "\tLikelihood = "
						+ String.format("%.8f", newLikelihood.getLogValue()) + "\tRates:" + ratesStr.toString());
				maxLikelihood = newLikelihood;
			}

		}
		if (errorLog != null)
			errorLog.write(logStr.toString());
	}

	// Rate Parameters Estimation
	// Maximum likelihood estimation available at
	// http://en.wikipedia.org/wiki/Gamma_distribution
	public static Pair<Double, Double> getMLERates(double[] rates) {

		double x;
		double sum = 0.0;
		double logSum = 0.0;
		int N = rates.length;
		for (int idx = 0; idx < N; idx++) {

			if (rates[idx] == 0.0)
				return null;

			x = rates[idx];
			sum += x;
			logSum += Math.log(x);
		}

		double s = Math.log(sum) - Math.log(N) - (1.0 / N) * logSum;
		double k = (3.0 - s + Math.sqrt((s - 3.0) * (s - 3.0) + 24 * s)) / (12.0 * s);
		double theta = (1.0 / (k * N)) * sum;

		double mean = k * theta;
		double cv = Math.sqrt(k) / k;

		Pair<Double, Double> mleRates = new Pair<Double, Double>(mean, cv);

		return mleRates;
	}

	// This method will derive Heuristic based combine msa for genes to construct
	// the gene tree
	public static MSAData getGeneMSA(MSAData[] msa, Mapping[] domain2gene, Mapping gene2species, SequenceType seqType,
			PRNG prng) {

		int numberOfDomains = domain2gene.length;
		LinkedHashMap<String, String> seqMap = new LinkedHashMap<String, String>();
		HashMap<String, ArrayList<String>> domainsOnGene = new HashMap<String, ArrayList<String>>();

		// Initialize the maps
		for (String key : gene2species.map.keySet()) {
			domainsOnGene.put(key, new ArrayList<String>());
			seqMap.put(key, new String());
		}

		// For each domain
		for (int r = 0; r < numberOfDomains; r++) {

			// Existence the rth domain(s) mapped to genes
			for (Map.Entry<String, String> entry : domain2gene[r].map.entrySet()) {
				String domainID = entry.getKey();
				String geneID = entry.getValue();
				domainsOnGene.get(geneID).add(domainID);
			}

			// Concatinating the rth domain on that gene
			for (String geneID : domainsOnGene.keySet()) {

				ArrayList<String> domainID = domainsOnGene.get(geneID);
				StringBuilder seq = new StringBuilder();

				/*
				 * Three possibilities: 1. Indel of rth domain, hence padding with X to maintain
				 * the alignment size) 2. Single instance of rth domain to concatinate 3. More
				 * than one domains, therefore picked randomly
				 */

				if (domainID.size() == 0) {
					while (seq.length() < msa[r].getNoOfPositions())
						seq.append("X");
				} else if (domainID.size() == 1) {
					seq.append(msa[r].getSequence(domainID.get(0)));
				} else {
					int randomChoice = prng.nextInt(domainID.size());
					seq.append(msa[r].getSequence(domainID.get(randomChoice)));

					// int fixedChoice = 0 ;
					// if(r==1) fixedChoice = 4;
					// seq.append( msa[r].getSequence( domainID.get( fixedChoice )));
				}

				// Now saving sequence as an alignment.
				String prevSeq = seqMap.get(geneID);
				seqMap.put(geneID, prevSeq + seq.toString());

			}

			// Clearing for next iteration.
			for (String key : domainsOnGene.keySet())
				domainsOnGene.get(key).clear();

		}

		// Creating MSAData object to return.
		Sequences geneAlignment = new Sequences();
		for (String identifier : seqMap.keySet())
			geneAlignment.append(identifier, seqMap.get(identifier));
		MSAData geneMSA = new MSAData(seqType, geneAlignment);

		return geneMSA;
	}

}

////////////////////////////////////////////////////
/*
 * 
 * // This method will initialize the paramters with actual simulation
 * parameters. public static void InitializeRatesWithRandomValues(Parameters pi)
 * {
 * 
 * // Birth death parameters pi.geneBDRates[0] = Math.random();
 * pi.geneBDRates[1] = Math.random(); for(int i=0 ; i<pi.numberOfDomains ; i++)
 * { pi.domBDRates[i][0] = Math.random(); pi.domBDRates[i][1] = Math.random(); }
 * // Edge rate parameters for(int i=0 ; i<pi.numberOfDomains ; i++) {
 * pi.domEdgeRates[i][0] = 4.0; pi.domEdgeRates[i][1] = 2.0; } }
 * 
 * ///////////////////////////////////////////////////// public static MSAData
 * getGeneAlignment(MSAData [] msa,Mapping [] domain2gene, Mapping
 * gene2species,SubstitutionMatrix Q, String geneAlignmentFile) throws
 * IOException {
 * 
 * 
 * LinkedHashMap<String,String> seqMap = new LinkedHashMap<String,String>();
 * HashMap<String, ArrayList<String>> domainsOnGene = new HashMap<String,
 * ArrayList<String>>();
 * 
 * // Initialize the maps for(String key: gene2species.map.keySet() ) {
 * domainsOnGene.put(key, new ArrayList<String>()); seqMap.put(key, new
 * String()); }
 * 
 * // For each domain for(int i=0 ; i< domain2gene.length ; i++ ) {
 * 
 * // Initialize the distance matrices double[][] distance =
 * DistanceMatrix.computeJC69(msa[i]);
 * 
 * // Determine the domain(s) mapped to a gene for (Map.Entry<String, String>
 * entry : domain2gene[i].map.entrySet()) { String domainID = entry.getKey();
 * String geneID = entry.getValue(); domainsOnGene.get(geneID).add(domainID); }
 * 
 * 
 * // Concatinating the domains on that gene for(String geneID :
 * domainsOnGene.keySet() ) {
 * 
 * ArrayList<String> domains = domainsOnGene.get(geneID); StringBuilder seq =
 * new StringBuilder();
 * 
 * if( domains.size() == 0 ) { // Padding with X to maintain the alignment size
 * while( seq.length() <= msa[i].getNoOfPositions() ) seq.append("X"); } else if
 * (domains.size() == 1 ) { // Single domain for concatination String idx
 * =domains.get(0); seq.append( msa[i].getSequence(idx) ); } else { // Pick the
 * best representative
 * 
 * double sumOfDistances [] = new double[domains.size()]; ArrayList <Integer>
 * domSet = new ArrayList <Integer>();
 * 
 * // Initialize the domain set for( int ii=0 ; ii < domains.size() ; ii++ ) {
 * Integer idx = (Integer) msa[i].getSequenceIndex(domains.get(ii));
 * domSet.add(idx); }
 * 
 * // Find the distances to other domain for( int ii=0 ; ii< domains.size() ;
 * ii++ ) { double sum=0 ; int idx = msa[i].getSequenceIndex(domains.get(ii));
 * for(int jj=0 ; jj< distance.length ; jj++) { if(!domSet.contains(jj)) sum +=
 * distance[idx][jj]; } sumOfDistances[ii] = sum; }
 * 
 * // Picking best one based on minimum distance criteria. int bestOne =
 * Mathematics.minIndex(sumOfDistances); seq.append( msa[i].getSequence(
 * domains.get( bestOne ) ) ); //System.out.println( geneID + "\t" +
 * domains.get( bestOne )); }
 * 
 * // Now saving sequence as an alignment. String prevSeq = seqMap.get(geneID);
 * seqMap.put(geneID, prevSeq + seq.toString() ); }
 * 
 * // Clearing for next iteration. for(String key: domainsOnGene.keySet() )
 * domainsOnGene.put(key, new ArrayList<String>()) ;
 * 
 * }
 * 
 * // Creating MSAData object to return. Sequences geneAlignment = new
 * Sequences(); for(String identifier: seqMap.keySet())
 * geneAlignment.append(identifier, seqMap.get(identifier));
 * 
 * // Writing result to file for validation purpose if(geneAlignmentFile!=null)
 * FastaRW.write(geneAlignmentFile, geneAlignment);
 * 
 * MSAData geneMSA = new MSAData(Q.getSequenceType(),geneAlignment);
 * 
 * 
 * 
 * return geneMSA; }
 */

///////////////////////////////
// double time = pi.domainTree[r].at[idx];
// geneDLModelLikelihood =
/////////////////////////////// TopologyPriors.getBirthDeathTreePriorFull(geneRecon,pi.geneBDRates,speciesExtinctionProb,null);
// rateModelLikelihoods =
/////////////////////////////// RateModel.getLogLikelihood1(pi.domEdgeRates[r],pi.domainTree[r],log);

/*
 * if( length<= 0 && time == 0.0 ) { logStr.
 * append("\nError: Either Branch length less than or equal to zero or arc time is zero"
 * ); continue; } else
 */

// int pIdx = pi.domainTree[r].nodes.get(idx).parent.id;
// double time = pi.domainTree[r].vt[pIdx] - pi.domainTree[r].vt[idx];
/*
 * import seqEvolution.LikelihoodVectors; public static void
 * estimateLengths(Parameters pi) {
 * 
 * int iteration=0; int maxIteration = 10;
 * 
 * Tree tree = pi.domainTree[0];
 * 
 * SubstitutionModel smodel;
 * 
 * LogDouble currentDataLikelihood = new LogDouble(1.0); LogDouble
 * oldDataLikelihood = new LogDouble(0.0);
 * 
 * do { smodel = new SubstitutionModel(pi.domMSA[0], pi.Q, tree,false);
 * 
 * // Compute likelihood of the tree currentDataLikelihood =
 * smodel.getLogLikelihood(); System.out.println("Current Data likelihood = " +
 * currentDataLikelihood.getLogValue());
 * 
 * 
 * //Estimating MLE ancestral sequences LinkedHashMap<Integer, String>
 * mleAncestralSeqs = smodel.getMLEAncestralSeq();
 * 
 * 
 * //Based on these sequences compute mle lengths ArrayList<Node> preOrder = new
 * ArrayList<Node>(); tree.getTreePreOrder(preOrder, tree.root);
 * 
 * String seq1,seq2; double [] lengths = new double[tree.nnodes];
 * 
 * for(Node node: preOrder) { for(int i=0 ; i<node.nchildren ; i++) { Node child
 * = node.children.get(i); seq1 = mleAncestralSeqs.get(node.id); seq2 =
 * mleAncestralSeqs.get(child.id); lengths[child.id] =
 * DistanceMatrix.computeMLE_BR(seq1,seq2, pi.Q); } }
 * 
 * // Now setting new lengths tree.setLengths(lengths);
 * 
 * 
 * }while(currentDataLikelihood.greaterThanOrEquals(oldDataLikelihood) &&
 * iteration++ < maxIteration);
 * 
 * }
 * 
 * public LinkedHashMap<Integer, String> getMLEAncestralSeq() {
 * 
 * LinkedHashMap<Integer, String> mleAncestralSeq = new LinkedHashMap<Integer,
 * String>();
 * 
 * this.computeModelLikelihood();
 * 
 * for(int i=0 ; i<this.T.nnodes ; i++) {
 * 
 * Node node = this.T.nodes.get(i); if (node.isLeaf()) {
 * mleAncestralSeq.put(node.id, this.D.getSequence(node.name)); } else { double
 * maxlikelihood; int MaxAlphabet; StringBuilder mleSeq = new StringBuilder();
 * LikelihoodVectors pl = this.mylikelihoods.get(node.id); for(int m=0 ; m <
 * this.D.getNoOfPositions() ; m++) { maxlikelihood = Double.NEGATIVE_INFINITY;
 * MaxAlphabet = 0; for(int a=0; a<this.alphabetSize ; a++) { if(
 * pl.likelihoods[m].get(a) > maxlikelihood ) { maxlikelihood =
 * pl.likelihoods[m].get(a); MaxAlphabet = a; } }
 * mleSeq.append(this.Q.getSequenceType().int2char(MaxAlphabet)); }
 * assert(mleSeq.length()!=this.D.getNoOfPositions());
 * mleAncestralSeq.put(node.id, mleSeq.toString()); } }
 * 
 * return mleAncestralSeq; }
 * 
 */
