package Main;

import java.io.IOException;
import java.math.BigInteger;

import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.analysis.solvers.UnivariateSolver;
import org.ejml.alg.dense.decomposition.DecompositionFactory;
import org.ejml.alg.dense.decomposition.EigenDecomposition;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import common.Log;
import common.LogDouble;
import common.PRNG;

import dlModel.BirthDeathProbs;
import dlModel.DLModel;

import phylogeny.Mapping;
import phylogeny.NewickReader;
import phylogeny.Reconciliation;
import phylogeny.Tree;
import rateModel.GammaDistribution;
import rateModel.RateModel;
import seqEvolution.AdditionalEJMLOps;
import seqEvolution.FastaRW;
import seqEvolution.MSAData;
import seqEvolution.RKSeqEvoModel;
import seqEvolution.SMD;
import seqEvolution.Sequences;
import seqEvolution.SubstitutionMatrix;
import seqEvolution.SubstitutionModel;
import treeConstructor.DistanceMatrix;
import treeConstructor.LikelihoodFunction;

public class testSystem {

	public static String CWD = "/Users/auwnm/Documents/Jworkspace/JPrIME_runs/GenPhyloData/Simulator/sim5/";
	public static Log log;
	public static PRNG prng;

	public static void main(String[] args) {

		try {

			System.out.println("Hello");

			String seed = null;
			log = new Log(CWD + "system.inf");
			prng = new PRNG(new BigInteger(seed));

			// actualDataLikelihood();
			// testRKSeqEvolutionModel() ;
			// testSeqEvolutionModel();

		} catch (Exception e) {
			System.out.println(e.toString());
			e.printStackTrace();
		}
	}

	public static void actualDataLikelihood() throws IOException {

		LogDouble overall = new LogDouble(1.0);
		overall.mult(testDLModel());
		overall.mult(testRTModel());
		overall.mult(testSeqEvolutionModel());

		System.out.println("\nOverall likelihood by all models : " + overall.getLogValue());
		System.out.println("=================================================================");

	}

	public static LogDouble testDLModel() throws IOException {

		int numberOfDomains = 2;
		String geneMapFile;
		String[] domMappingFiles = new String[numberOfDomains];
		String[] domTreeFile = new String[numberOfDomains];

		boolean useRootEdge = true;

		String speciesTreeFile = CWD + "hostTree.nwk";
		String geneTreeFile = CWD + "geneTree.nwk";
		domTreeFile[0] = CWD + "domainTree1.nwk";
		domTreeFile[1] = CWD + "domainTree2.nwk";
		geneMapFile = CWD + "geneTree.map";
		domMappingFiles[0] = CWD + "domainTree1.map";
		domMappingFiles[1] = CWD + "domainTree2.map";

		// Reading species tree
		Tree speciesTree = NewickReader.readNewickTreeFile(speciesTreeFile);
		speciesTree.setName("speciesTree");
		speciesTree.setDepths(speciesTree.root, 0);
		speciesTree.setTimesFromOriginalLengths();

		// Reading gene tree
		Tree geneTree = NewickReader.readNewickTreeFile(geneTreeFile);
		geneTree.setName("geneTree");
		geneTree.setDepths(geneTree.root, 0);
		geneTree.setTimesFromOriginalLengths();

		// Reading domain trees
		Tree domainTree[] = new Tree[numberOfDomains];
		domainTree[0] = NewickReader.readNewickTreeFile(domTreeFile[0]);
		domainTree[0].setName("domainTree1");
		domainTree[1] = NewickReader.readNewickTreeFile(domTreeFile[1]);
		domainTree[1].setName("domainTree2");

		// Reading mapping files
		Mapping gene2species = new Mapping();
		Mapping[] domain2genes = new Mapping[2];
		gene2species.readMapping(geneMapFile);
		domain2genes[0] = new Mapping();
		domain2genes[0].readMapping(domMappingFiles[0]);
		domain2genes[1] = new Mapping();
		domain2genes[1].readMapping(domMappingFiles[1]);

		// Computing MPR reconciliations between domain,gene and species tree
		Reconciliation geneRecon = new Reconciliation(geneTree, speciesTree, gene2species);
		Reconciliation[] domainRecon = new Reconciliation[numberOfDomains];
		for (int i = 0; i < numberOfDomains; i++)
			domainRecon[i] = new Reconciliation(domainTree[i], geneTree, domain2genes[i]);

		// Specfying birth and death rate for gene, domain trees ...
		double[] geneBDRates = { 0.012, 0.01 };
		double[] dom1BDRates = { 0.012003615764595185, 0.010260696634929381 };
		double[] dom2BDRates = { 0.012700135230749, 0.010995735881647938 };

		// Computing Extinction Prob
		double[] speciesExtinctionProb = BirthDeathProbs.calcExtinctionProb(speciesTree, geneBDRates);
		double[] gene1ExtinctionProb = BirthDeathProbs.calcExtinctionProb(geneTree, dom1BDRates);
		double[] gene2ExtinctionProb = BirthDeathProbs.calcExtinctionProb(geneTree, dom2BDRates);

		DLModel geneDLModel = new DLModel(geneRecon, geneBDRates, speciesExtinctionProb, useRootEdge, prng, log);
		DLModel dom1DLModel = new DLModel(domainRecon[0], dom1BDRates, gene1ExtinctionProb, useRootEdge, prng, log);
		DLModel dom2DLModel = new DLModel(domainRecon[1], dom2BDRates, gene2ExtinctionProb, useRootEdge, prng, log);

		// Likelihood computation
		LogDouble geneDLModelLikelihood = geneDLModel.birthDeathTreePrior();// TopologyPriors.getBirthDeathTreePrior(geneRecon,geneBDRates,speciesExtinctionProb,null);
		LogDouble dom1DLModelLikelihood = dom1DLModel.birthDeathTreePrior();// TopologyPriors.getBirthDeathTreePrior(domainRecon[0],dom1BDRates,gene1ExtinctionProb,null);
		LogDouble dom2DLModelLikelihood = dom2DLModel.birthDeathTreePrior();// TopologyPriors.getBirthDeathTreePrior(domainRecon[1],dom2BDRates,gene2ExtinctionProb,null);

		// Computing overall likelihood
		LogDouble overall = new LogDouble(1.0);
		overall.mult(geneDLModelLikelihood);
		overall.mult(dom1DLModelLikelihood);
		overall.mult(dom2DLModelLikelihood);

		// Print likelihoods for each tree
		System.out.println("\nDL Model:");
		System.out.println(
				"Gene tree     likelihood P( T_G , r | lambda,mu,T_S ) = " + geneDLModelLikelihood.getLogValue());
		System.out.println(
				"Domain tree 1 likelihood P( T_D1, r | lambda_D1,mu_D2,T_G ) = " + dom1DLModelLikelihood.getLogValue());
		System.out.println(
				"Domain tree 1 likelihood P( T_D2, r | lambda_D1,mu_D2,T_G ) = " + dom2DLModelLikelihood.getLogValue());
		System.out.println("Overall likelihood by DL Model = " + overall.getLogValue());

		return overall;

	}

	public static LogDouble testRTModel() throws IOException {

		int numberOfDomains = 2;
		String[] domainTreeLengthFile = new String[numberOfDomains];
		String[] domainTreeTimeFile = new String[numberOfDomains];

		boolean useRootEdge = true;

		domainTreeLengthFile[0] = CWD + "relaxed_domainTree1.nwk";
		domainTreeLengthFile[1] = CWD + "relaxed_domainTree2.nwk";
		domainTreeTimeFile[0] = CWD + "domainTree1.nwk";
		domainTreeTimeFile[1] = CWD + "domainTree2.nwk";

		// Reading original domain trees with lengths
		Tree domainTreeWithLengths[] = new Tree[2];
		domainTreeWithLengths[0] = NewickReader.readNewickTreeFile(domainTreeLengthFile[0]);
		domainTreeWithLengths[0].setName("domainTree1");
		domainTreeWithLengths[0].setBranchLengthsFromOriginalLengths();
		domainTreeWithLengths[1] = NewickReader.readNewickTreeFile(domainTreeLengthFile[1]);
		domainTreeWithLengths[1].setName("domainTree2");
		domainTreeWithLengths[1].setBranchLengthsFromOriginalLengths();

		// Reading original domain trees with time information
		Tree domainTree[] = new Tree[numberOfDomains];
		for (int i = 0; i < numberOfDomains; i++) {
			domainTree[i] = NewickReader.readNewickTreeFile(domainTreeTimeFile[i]);
			domainTree[i].setName("domainTree" + (i + 1));
			domainTree[i].setTimesFromOriginalLengths();

			// Setting time info from originals
			domainTreeWithLengths[i].vt = domainTree[i].getNormalizedTimes();
			domainTreeWithLengths[i].peakTime = 1.0;
		}

		// Edge rate parameters
		double k = 0.49;
		double theta = 1.1;
		double[] dom1EdgeRates = { k * theta, Math.sqrt(k) / k };
		double[] dom2EdgeRates = { k * theta, Math.sqrt(k) / k };

		// Calling Rate model ....
		// Rate Model initialized
		RateModel rateModel1 = new RateModel(dom1EdgeRates, useRootEdge); // domainTreeWithLengths[0]
		RateModel rateModel2 = new RateModel(dom2EdgeRates, useRootEdge); // domainTreeWithLengths[1]

		rateModel1.update(domainTreeWithLengths[0]);
		rateModel2.update(domainTreeWithLengths[1]);

		LogDouble rateModelLikelihoods[] = new LogDouble[numberOfDomains];
		rateModelLikelihoods[0] = rateModel1.getLogLikelihood();
		rateModelLikelihoods[1] = rateModel2.getLogLikelihood();

		// Computing overall likelihood
		LogDouble overall = new LogDouble(1.0);
		overall.mult(rateModelLikelihoods[0]);
		overall.mult(rateModelLikelihoods[1]);

		// Print likelihoods for each tree
		System.out.println("\nRate Model:");
		System.out.println("Rate Model likelihood for domain tree 1 : " + rateModelLikelihoods[0].getLogValue());
		System.out.println("Rate Model likelihood for domain tree 2 : " + rateModelLikelihoods[1].getLogValue());
		System.out.println("Overall likelihood by Rate Model : " + overall.getLogValue());

		return overall;

	}

	public static LogDouble testSeqEvolutionModel() throws IOException {

		// Specify the substitution model
		String substitutionModel = "JC69";
		SubstitutionMatrix Q = SMD.create(substitutionModel);

		// Reading original domain trees
		String dom1TreeFile = CWD + "relaxed_domainTree1.nwk";
		String dom2TreeFile = CWD + "relaxed_domainTree2.nwk";

		Tree domainTree[] = new Tree[2];
		domainTree[0] = NewickReader.readNewickTreeFile(dom1TreeFile);
		domainTree[0].setName("domainTree1");
		domainTree[0].setBranchLengthsFromOriginalLengths();

		domainTree[1] = NewickReader.readNewickTreeFile(dom2TreeFile);
		domainTree[1].setName("domainTree2");
		domainTree[1].setBranchLengthsFromOriginalLengths();

		// Reading MSAs
		String[] alignmentFiles = new String[2];
		alignmentFiles[0] = CWD + "seq_relaxed_domainTree1.fasta";
		alignmentFiles[1] = CWD + "seq_relaxed_domainTree2.fasta";
		MSAData[] domMSA = readAlignmentFiles(alignmentFiles, Q);

		// Computing Likelihood
		SubstitutionModel smodel[] = new SubstitutionModel[2];
		smodel[0] = new SubstitutionModel(domMSA[0], Q, domainTree[0], true);
		smodel[1] = new SubstitutionModel(domMSA[1], Q, domainTree[1], true);

		// Computing overall likelihood
		LogDouble overall = new LogDouble(1.0);

		// To estimate the performance
		long startTime = System.nanoTime();

		overall.mult(smodel[0].getLogLikelihood());
		overall.mult(smodel[1].getLogLikelihood());

		long stopTime = System.nanoTime();
		System.out.println("\nApprox. computational time (sec) = " + (stopTime - startTime) / Math.pow(10, 9));

		// Print likelihoods for each tree
		System.out.println("\nSequence Evolution Model:");
		System.out.println("Data likelihood for domain tree 1 : " + smodel[0].getLogLikelihood().getLogValue());
		System.out.println("Data likelihood for domain tree 2 : " + smodel[1].getLogLikelihood().getLogValue());
		System.out.println("Overall data likelihood for domain trees :  " + overall.getLogValue());

		return overall;
	}

	public static LogDouble testRKSeqEvolutionModel() throws IOException {

		// Loading related libaray
		System.load("/Users/auwnm/Documents/cworkspace/RasmussenSeqEvo/Debug/libRasmussenSeqEvo.dylib");

		// Specify the substitution model
		String substitutionModel = "JC69";
		SubstitutionMatrix Q = SMD.create(substitutionModel);

		// Reading original domain trees
		String dom1TreeFile = CWD + "relaxed_domainTree1.nwk";
		String dom2TreeFile = CWD + "relaxed_domainTree2.nwk";

		Tree domainTree[] = new Tree[2];
		domainTree[0] = NewickReader.readNewickTreeFile(dom1TreeFile);
		domainTree[0].setName("domainTree1");
		domainTree[0].setBranchLengthsFromOriginalLengths();

		domainTree[1] = NewickReader.readNewickTreeFile(dom2TreeFile);
		domainTree[1].setName("domainTree2");
		domainTree[1].setBranchLengthsFromOriginalLengths();

		// Reading MSAs
		String[] alignmentFiles = new String[2];
		alignmentFiles[0] = CWD + "seq_relaxed_domainTree1.fasta";
		alignmentFiles[1] = CWD + "seq_relaxed_domainTree2.fasta";
		MSAData[] domMSA = readAlignmentFiles(alignmentFiles, Q);

		// To estimate the performance
		long startTime = System.nanoTime();

		LogDouble likelihood1 = RKSeqEvoModel.getSeqLikelihood(domainTree[0], domMSA[0]);
		LogDouble likelihood2 = RKSeqEvoModel.getSeqLikelihood(domainTree[1], domMSA[1]);

		long stopTime = System.nanoTime();
		System.out.println("\nApprox. computational time (sec) = " + (stopTime - startTime) / Math.pow(10, 9));

		// Print likelihoods for each tree
		System.out.println("\nSequence Evolution Model:");
		System.out.println("Data likelihood for domain tree 1 : " + likelihood1.getLogValue());
		System.out.println("Data likelihood for domain tree 2 : " + likelihood2.getLogValue());

		return null;
	}

	// Reading alignment files
	public static MSAData[] readAlignmentFiles(String[] alignmentFiles, SubstitutionMatrix Q) throws IOException {

		// Reading domain alignment files one by one
		int noOfAlignmentFiles = alignmentFiles.length;
		Sequences[] alignments = new Sequences[noOfAlignmentFiles];
		for (int i = 0; i < noOfAlignmentFiles; i++) {
			alignments[i] = FastaRW.readFastaAlign(alignmentFiles[i]);
			if (alignments[i] == null) {
				return null;
			}
		}
		// Making MSA objects from these alignments
		MSAData[] domainMSA = new MSAData[noOfAlignmentFiles];
		for (int i = 0; i < noOfAlignmentFiles; i++)
			domainMSA[i] = new MSAData(Q.getSequenceType(), alignments[i]);

		return domainMSA;
	}

	// For debugging purposes ...

	// Imp Optimizing hill climbing ....

	/*
	 * 
	 * 
	 * package MCMC;
	 * 
	 * import java.io.IOException; import java.util.ArrayList; import
	 * java.util.LinkedHashMap;
	 * 
	 * import phylogeny.Node; import phylogeny.Reconciliation; import
	 * rateModel.RateModel; import seqEvolution.SubstitutionModel; import
	 * Main.Parameters; import Main.ParametersHandler; import Main.mainSystem;
	 * 
	 * import common.Log; import common.LogDouble; import dlModel.BirthDeathProbs;
	 * import dlModel.TopologyPriors;
	 * 
	 * public class HillClimbingOnState {
	 * 
	 * static double RESOLUTION = 1e5; static double ACCELERATION = 1.2; static
	 * double MAXTRIES = 1; static boolean WriteStateParameters = false;
	 * 
	 * public static LogDouble compute(Parameters pi, Log log) throws IOException {
	 * 
	 * // To estimate the performance long startTime = System.nanoTime();
	 * 
	 * 
	 * // Computing MPR reconciliations between domain,gene and species tree
	 * Reconciliation geneRecon = new
	 * Reconciliation(pi.geneTree,pi.speciesTree,pi.gene2species); Reconciliation []
	 * domainRecon = new Reconciliation[pi.numberOfDomains]; for(int i=0 ; i <
	 * pi.numberOfDomains ; i++) domainRecon[i] = new
	 * Reconciliation(pi.domainTree[i],pi.geneTree, pi.domain2genes[i]);
	 * 
	 * 
	 * // Writing parameters if(WriteStateParameters)
	 * log.writeParameters(pi,geneRecon,domainRecon,true);
	 * 
	 * 
	 * LogDouble stateLikelihood = new LogDouble (1.0);
	 * 
	 * 
	 * double [] extinctionProb =
	 * BirthDeathProbs.calcExtinctionProb(pi.speciesTree,pi.geneBDRates); LogDouble
	 * geneDLModelLikelihood = TopologyPriors.getBirthDeathTreePrior(geneRecon,
	 * pi.geneBDRates, extinctionProb, log);
	 * stateLikelihood.mult(geneDLModelLikelihood);
	 * 
	 * 
	 * //Random-restart hill climbing LogDouble bestLikelihood,currentLikelihood;
	 * 
	 * int trial=0; bestLikelihood = new LogDouble(0.0); for( ; trial< MAXTRIES ;
	 * trial ++ ) {
	 * 
	 * // Initialize the arcs and vertices times on gene and domain tree
	 * TopologyPriors.sampleRealisation(geneRecon, pi.geneBDRates, extinctionProb,
	 * log); for(int r=0 ; r < pi.numberOfDomains ; r++) { extinctionProb =
	 * BirthDeathProbs.calcExtinctionProb(pi.geneTree,pi.domBDRates[r]);
	 * TopologyPriors.sampleRealisation(domainRecon[r], pi.domBDRates[r],
	 * extinctionProb,log); }
	 * 
	 * currentLikelihood = maximizeGeneLikelihood(pi ,geneRecon,domainRecon,log);
	 * 
	 * if(currentLikelihood.greaterThan(bestLikelihood)) bestLikelihood =
	 * currentLikelihood;
	 * 
	 * }
	 * 
	 * // Updating overall likelihood stateLikelihood.mult(bestLikelihood);
	 * 
	 * long stopTime = System.nanoTime(); System.out.println(
	 * "\nApprox. computational time (sec) = " + (double) (stopTime - startTime) /
	 * Math.pow(10,9));
	 * 
	 * return stateLikelihood; }
	 * 
	 * 
	 * 
	 * public static LogDouble maximizeGeneLikelihood(Parameters pi ,Reconciliation
	 * geneRecon,Reconciliation [] domainRecon,Log log) throws IOException {
	 * 
	 * int iteration=0; double delta; double stepSize; int bestStrategy; boolean
	 * moved,recovered; LogDouble bestLikelihood,currentLikelihood,tempLikelihood;
	 * double [] candidate = {-ACCELERATION , -1.0/ACCELERATION, 0.0,
	 * 1.0/ACCELERATION, ACCELERATION }; stepSize = geneRecon.host.getPeakTime() /
	 * RESOLUTION ;
	 * 
	 * 
	 * Log stateLog = new Log(mainSystem.CWD + "log.state.out");
	 * 
	 * 
	 * LogDouble initialLikelihood = new LogDouble(0.0);
	 * 
	 * // Outermost loop while(true) {
	 * 
	 * 
	 * HillState[] domainTreeTimes = new HillState[pi.numberOfDomains]; for(int r=0
	 * ; r<pi.numberOfDomains ; r++) domainTreeTimes[r] = new
	 * HillState(domainRecon[r]);
	 * 
	 * currentLikelihood = new LogDouble(0.0); ArrayList<Node> geneDuplications =
	 * geneRecon.getDuplications(); for(Node node: geneDuplications) { delta =
	 * stepSize * getBranchWeight(node,geneRecon,geneRecon.host.getLargerArcTime());
	 * do { // Likelihood of the current state currentLikelihood =
	 * maximizeDomainLikelihood(pi ,domainRecon,log);
	 * 
	 * // Saving current times for domain Trees for(int r=0 ; r<pi.numberOfDomains ;
	 * r++) domainTreeTimes[r].saveTimeState();
	 * 
	 * // Picking the best possible state bestStrategy = -1; bestLikelihood = new
	 * LogDouble(0.0); for(int c=0 ; c<5 ; c++) {
	 * 
	 * moved = moveGeneNode(node,delta * candidate[c],geneRecon, domainRecon);
	 * if(moved) { tempLikelihood = maximizeDomainLikelihood(pi ,domainRecon,log);
	 * if(tempLikelihood.greaterThanOrEquals(bestLikelihood)) { bestLikelihood =
	 * tempLikelihood; bestStrategy = c; } recovered = moveGeneNode(node, -1.0 *
	 * delta * candidate[c],geneRecon, domainRecon); assert(recovered); for(int r=0
	 * ; r<pi.numberOfDomains ; r++) domainTreeTimes[r].recoverTimeState(); } }
	 * 
	 * // Commit best move and accelerate if( candidate[bestStrategy] != 0.0 ) {
	 * moveGeneNode(node, delta * candidate[bestStrategy], geneRecon, domainRecon);
	 * delta *= candidate[bestStrategy]; stateLog.writeHillIteration(iteration++,
	 * node.id, delta * candidate[bestStrategy], bestLikelihood,currentLikelihood);
	 * }
	 * 
	 * }while(bestLikelihood.greaterThan(currentLikelihood));
	 * 
	 * }
	 * 
	 * 
	 * if(currentLikelihood.greaterThan(initialLikelihood)) initialLikelihood =
	 * currentLikelihood; else break; }
	 * 
	 * stateLog.close(); return currentLikelihood; }
	 * 
	 * 
	 * public static LogDouble maximizeDomainLikelihood(Parameters pi
	 * ,Reconciliation [] domainRecon,Log log) throws IOException {
	 * 
	 * ArrayList<Node> domDuplications; LogDouble domDLModelLikelihood; double []
	 * extinctionProb; LogDouble domainLikelihood = new LogDouble (1.0);
	 * 
	 * double delta; double stepSize; int bestStrategy; boolean moved,recovered;
	 * LogDouble bestLikelihood,currentLikelihood,tempLikelihood; double []
	 * candidate = {-ACCELERATION , -1.0/ACCELERATION, 0.0, 1.0/ACCELERATION,
	 * ACCELERATION };
	 * 
	 * 
	 * // Assuming time on gene tree has changed. for(int r=0 ; r <
	 * pi.numberOfDomains ; r++) {
	 * 
	 * // DL-Model Likelihood for this domain tree extinctionProb =
	 * BirthDeathProbs.calcExtinctionProb(pi.geneTree,pi.domBDRates[r]);
	 * domDLModelLikelihood = TopologyPriors.getBirthDeathTreePrior(domainRecon[r],
	 * pi.domBDRates[r], extinctionProb,log);
	 * domainLikelihood.mult(domDLModelLikelihood);
	 * 
	 * 
	 * LogDouble initialLikelihood = new LogDouble(0.0);
	 * 
	 * // Entering in outermost loop while(true) {
	 * 
	 * // Optimisation for each duplication node of this domain tree stepSize =
	 * domainRecon[r].host.getPeakTime() / RESOLUTION ; currentLikelihood = new
	 * LogDouble(0.0); domDuplications = domainRecon[r].getDuplications(); for(Node
	 * node: domDuplications) { delta = stepSize *
	 * getBranchWeight(node,domainRecon[r],domainRecon[r].host.getLargerArcTime());
	 * do { // Likelihood of the current state currentLikelihood =
	 * RateModel.getLogLikelihood(pi.domEdgeRates[r],pi.domainTree[r],log);
	 * 
	 * // Picking the best possible state bestStrategy = -1; bestLikelihood = new
	 * LogDouble(0.0); for(int c=0 ; c<5 ; c++) { moved = moveDomainNode(node, delta
	 * * candidate[c], domainRecon[r]); if(moved) { tempLikelihood =
	 * RateModel.getLogLikelihood(pi.domEdgeRates[r],pi.domainTree[r],log);
	 * if(tempLikelihood.greaterThanOrEquals(bestLikelihood)) { bestLikelihood =
	 * tempLikelihood; bestStrategy = c; } recovered = moveDomainNode(node, -1.0 *
	 * delta * candidate[c], domainRecon[r]); assert(recovered); } }
	 * 
	 * // Commit this best move and accelerate if( candidate[bestStrategy] != 0.0 )
	 * { moveDomainNode(node, delta * candidate[bestStrategy], domainRecon[r]);
	 * delta *= candidate[bestStrategy]; }
	 * 
	 * }while(bestLikelihood.greaterThan(currentLikelihood)); }
	 * if(currentLikelihood.greaterThan(initialLikelihood)) initialLikelihood =
	 * currentLikelihood; else break;
	 * 
	 * }
	 * 
	 * // Updating overall likelihood domainLikelihood.mult(initialLikelihood); }
	 * 
	 * 
	 * 
	 * 
	 * return domainLikelihood; }
	 * 
	 * 
	 * 
	 * private static boolean moveGeneNode(Node gnode,double delta,Reconciliation
	 * geneRecon,Reconciliation [] domainRecon ) {
	 * 
	 * int numberOfDomains = domainRecon.length;
	 * LinkedHashMap<Integer,ArrayList<Node>> influencedNodes = new
	 * LinkedHashMap<Integer,ArrayList<Node>>();
	 * 
	 * 
	 * double [] validBounds = getBounds(gnode,geneRecon);
	 * 
	 * for(int r=0 ; r < numberOfDomains ; r++ ) {
	 * 
	 * ArrayList<Node> influenced = new ArrayList<Node>(); for(int nid=0; nid <
	 * domainRecon[r].guestPrime.nnodes ; nid++) { Node dnode =
	 * domainRecon[r].guestPrime.nodes.get(nid); if(
	 * domainRecon[r].nodeMapPrime.get(dnode.id) == gnode.id &&
	 * domainRecon[r].eventMapPrime.get(dnode.id) == Reconciliation.EVENT_SPEC ) {
	 * 
	 * double [] bounds = getBounds(dnode,domainRecon[r]);
	 * 
	 * if(bounds[0] > validBounds[0]) validBounds[0] = bounds[0];
	 * 
	 * if(bounds[1] < validBounds[1]) validBounds[1] = bounds[1];
	 * 
	 * 
	 * if( dnode.id <= domainRecon[r].guest.root.id ) influenced.add(dnode); } }
	 * influencedNodes.put(r, influenced); }
	 * 
	 * 
	 * double t = geneRecon.guestPrime.vt[gnode.id];
	 * 
	 * if(delta==0) return true;
	 * 
	 * if( t+delta >= validBounds[1] || t+delta <= validBounds[0] ) return false;
	 * 
	 * 
	 * geneRecon.guestPrime.vt[gnode.id] = t+delta;
	 * 
	 * // Updating in original tree i.e tree without implied nodes // Assuming same
	 * duplication id in both the trees Node originalNode =
	 * geneRecon.guest.nodes.get(gnode.id); geneRecon.guest.vt[originalNode.id] =
	 * t+delta ; geneRecon.guest.at[originalNode.id] -= delta;
	 * geneRecon.guest.at[originalNode.children.get(0).id] += delta;
	 * geneRecon.guest.at[originalNode.children.get(1).id] += delta;
	 * 
	 * 
	 * 
	 * for(int r=0 ; r < numberOfDomains ; r++ ) { for(Node dnode:
	 * influencedNodes.get(r)) moveDomainNode(dnode,delta,domainRecon[r]); }
	 * 
	 * return true; }
	 * 
	 * private static boolean moveDomainNode(Node dnode,double delta,Reconciliation
	 * domainRecon) {
	 * 
	 * 
	 * double [] validBounds = getBounds(dnode,domainRecon);
	 * 
	 * 
	 * double t = domainRecon.guestPrime.vt[dnode.id];
	 * 
	 * if(delta==0) return true;
	 * 
	 * if( t+delta >= validBounds[1] || t+delta <= validBounds[0] ) return false;
	 * 
	 * 
	 * domainRecon.guestPrime.vt[dnode.id] = t+delta;
	 * 
	 * // Updating in original tree i.e tree without implied nodes // Assuming same
	 * duplication id in both the trees Node originalNode =
	 * domainRecon.guest.nodes.get(dnode.id); domainRecon.guest.vt[originalNode.id]
	 * = t+delta ; domainRecon.guest.at[originalNode.id] -= delta;
	 * domainRecon.guest.at[originalNode.children.get(0).id] += delta;
	 * domainRecon.guest.at[originalNode.children.get(1).id] += delta;
	 * 
	 * return true; }
	 * 
	 * private static double [] getBounds(Node node,Reconciliation recon) {
	 * 
	 * double [] bounds = new double [2];
	 * 
	 * double nodeTime = recon.guestPrime.vt[node.id];
	 * 
	 * double parentTime = node.isRoot() ? recon.host.getPeakTime() :
	 * recon.guestPrime.vt[node.parent.id] ;
	 * 
	 * double maxChildTime = 0.0; for(int i=0 ; i<node.nchildren ; i++) { double
	 * childTime = recon.guestPrime.vt[node.children.get(i).id]; if( childTime >
	 * maxChildTime ) maxChildTime = childTime; }
	 * 
	 * assert( nodeTime > maxChildTime && nodeTime < parentTime );
	 * 
	 * 
	 * bounds[0] = maxChildTime; bounds[1] = parentTime;
	 * 
	 * return bounds;
	 * 
	 * }
	 * 
	 * private static double getBranchWeight(Node node,Reconciliation R,double
	 * largerBranch) {
	 * 
	 * double t0;
	 * 
	 * if(!node.isRoot()) t0 = R.guestPrime.vt[node.parent.id]; else t0 =
	 * R.host.vt[R.host.root.id] + R.host.at[R.host.root.id];
	 * 
	 * double t1 = R.guestPrime.vt[node.children.get(0).id]; double t2 =
	 * R.guestPrime.vt[node.children.get(1).id]; double interval = t0 - Math.max(t1,
	 * t2);
	 * 
	 * assert( interval > Math.max(t1, t2) );
	 * 
	 * return interval/largerBranch; }
	 * 
	 * }
	 * 
	 * 
	 * 
	 * 
	 * 
	 */

	/*
	 * 
	 * public static LogDouble compute(Parameters pi, Log log) throws IOException {
	 * 
	 * 
	 * // Computing MPR reconciliations between domain,gene and species tree
	 * Reconciliation geneRecon = new
	 * Reconciliation(pi.geneTree,pi.speciesTree,pi.gene2species); Reconciliation []
	 * domainRecon = new Reconciliation[pi.numberOfDomains]; for(int i=0 ; i <
	 * pi.numberOfDomains ; i++) domainRecon[i] = new
	 * Reconciliation(pi.domainTree[i],pi.geneTree, pi.domain2genes[i]);
	 * 
	 * 
	 * // If this is the first MCMC iteration, initialize the birthdeath & rate
	 * parameters. if(FirstMCMCIteration) { ParametersHandler.estimateBDRates(pi,
	 * geneRecon, domainRecon,log);
	 * ParametersHandler.estimateEdgeRates(pi,geneRecon,domainRecon,log);
	 * //ParametersHandler.InitializeParametersWithActualEstimates(pi);
	 * log.writeRateParameters(pi);
	 * 
	 * }
	 * 
	 * // To estimate the performance long startTime = System.nanoTime();
	 * 
	 * 
	 * // Writing parameters //if(WriteStateParameters) //
	 * log.writeParameters(pi,geneRecon,domainRecon,true);
	 * 
	 * 
	 * LogDouble overallLikelihood = new LogDouble (1.0); double [] extinctionProb =
	 * BirthDeathProbs.calcExtinctionProb(pi.speciesTree,pi.geneBDRates); LogDouble
	 * geneDLModelLikelihood = TopologyPriors.getBirthDeathTreePrior(geneRecon,
	 * pi.geneBDRates, extinctionProb, log);
	 * overallLikelihood.mult(geneDLModelLikelihood);
	 * 
	 * 
	 * //Random-restart hill climbing LogDouble bestLikelihood,currentLikelihood;
	 * 
	 * int trial=0; bestLikelihood = new LogDouble(0.0); for( ; trial< MAXTRIES ;
	 * trial ++ ) {
	 * 
	 * // Initialize the arcs and vertices times on gene and domain tree
	 * TopologyPriors.sampleRealisation(geneRecon, pi.geneBDRates, extinctionProb,
	 * log); for(int r=0 ; r < pi.numberOfDomains ; r++) { extinctionProb =
	 * BirthDeathProbs.calcExtinctionProb(pi.geneTree,pi.domBDRates[r]);
	 * TopologyPriors.sampleRealisation(domainRecon[r], pi.domBDRates[r],
	 * extinctionProb,log); }
	 * 
	 * currentLikelihood = maximizeGeneLikelihood(pi ,geneRecon,domainRecon,log);
	 * 
	 * if(currentLikelihood.greaterThan(bestLikelihood)) bestLikelihood =
	 * currentLikelihood;
	 * 
	 * }
	 * 
	 * // Updating overall likelihood overallLikelihood.mult(bestLikelihood);
	 * System.out.println("Current Liklihood = " + overallLikelihood.getLogValue());
	 * 
	 * 
	 * // Performing sequence evolution at domain Level LogDouble dataLikelihood =
	 * new LogDouble(1.0); SubstitutionModel smodel[] = new
	 * SubstitutionModel[pi.numberOfDomains]; for(int i=0; i < pi.numberOfDomains ;
	 * i++ ) { smodel[i] = new SubstitutionModel(pi.domMSA[i], pi.Q,
	 * pi.domainTree[i],false); dataLikelihood.mult(smodel[i].getLogLikelihood());
	 * System.out.println("Data likelihood for domain tree " + (i+1) + ":" +
	 * smodel[i].getLogLikelihood()); }
	 * 
	 * 
	 * LogDouble stateLikelihood = new LogDouble(1.0);
	 * stateLikelihood.mult(dataLikelihood);
	 * stateLikelihood.mult(overallLikelihood);
	 * 
	 * 
	 * long stopTime = System.nanoTime(); System.out.println(
	 * "\nApprox. computational time (sec) = " + (double) (stopTime - startTime) /
	 * Math.pow(10,9));
	 * 
	 * return stateLikelihood; }
	 * 
	 * 
	 * 
	 * 
	 */

	/*
	 * 
	 * Length Proposer ...
	 * 
	 * // This method will perturb branch lengths w.r.t. its tree structure //
	 * nextInt = Returns a pseudorandom, uniformly distributed int value between 0
	 * (inclusive) and the specified value (exclusive) public boolean
	 * cacheAndPerturb() {
	 * 
	 * // First cache the current state of length parameter(s) this.lengthsCache =
	 * Arrays.copyOf(this.lengths, this.lengths.length);
	 * 
	 * 
	 * // Step 1: Randomly pick a path on a tree ArrayList<Node> path = new
	 * ArrayList<Node>(); pickPathRandomly(this.tree, null, path);
	 * 
	 * 
	 * // Step 2: Randomly pick branches along this path for perturbation int
	 * noOfPerturbed = this.prng.nextInt( path.size() ); for(int i=0;
	 * i<noOfPerturbed ; i++) { int picked = this.prng.nextInt( path.size() );
	 * perturbIdx.add(path.get(picked).id); }
	 * 
	 * // Initializing forward and backward this.forwardDensity = new
	 * LogDouble(1.0); this.backwardDensity = new LogDouble(1.0);
	 * 
	 * // Proposal is ready now Pair<double [],LogDouble[]> perturbedValues =
	 * NormalProposer.perturbMultiValues(this.lengths, perturbIdx, interval,
	 * this.proposalCV, prng);
	 * 
	 * 
	 * if(perturbedValues == null) return false;
	 * 
	 * // Changing the state of parameter for(int i=0 ; i<perturbIdx.size() ; i++)
	 * this.lengths[ perturbIdx.get(i) ] = perturbedValues.first[i];
	 * 
	 * this.forwardDensity = perturbedValues.second[0]; this.backwardDensity =
	 * perturbedValues.second[1];
	 * 
	 * return true; }
	 * 
	 * // Helping function for picking a path on a binary tree randomly public void
	 * pickPathRandomly(Tree tree, Node node, ArrayList<Node> path) { if(node ==
	 * null) node = tree.root;
	 * 
	 * path.add(node); if(node.isLeaf()) { return; } else { int childIdx =
	 * this.prng.nextInt(2);
	 * pickPathRandomly(tree,node.children.get(childIdx),path); } }
	 * 
	 * 
	 * 
	 * 
	 * 
	 */

	/*
	 * double x = 5.0; PRNG prng = new PRNG(); Pair<Double,LogDouble[]> rv =
	 * NormalProposer.perturbSingleValue(x,0.75,prng); System.out.println(rv.first +
	 * "  forward density = " + rv.second[0]+" Bacward Density =" +rv.second[1]);
	 * 
	 */
	/*
	 * // Logging gene tree sampled times log.writeTree(pi.geneTree, false, null,
	 * null);
	 * 
	 * // Logging extinction prob for domain tress log.write("\nFor domain " + (i+1)
	 * + " :"); log.writeExtinctionProb(pi.geneTree, extinctionProb);
	 * 
	 * // Logging domain tree sampled times log.writeTree(pi.domainTree[i], false,
	 * null, null);
	 */

	/*
	 * Parameters params = new
	 * Parameters(birthRate,deathRate,edgeRatePDMean,edgeRatePDCV);
	 * 
	 * System.out.println("Reading species tree"); Tree sTree =
	 * NewickReader.readNewickTreeFile(sTreeFile); sTree.setDepths(sTree.root,0);
	 * sTree.setVertexTimesFromArcTimes(); sTree.setNormalizeTimes();
	 * sTree.setSpeciationEvents(); sTree.show();
	 * 
	 * System.out.println("Computing Extinction Probs along host tree.");
	 * BirthDeathProbs.calcExtinctionProb(sTree,params);
	 * BirthDeathProbs.show(sTree);
	 * 
	 * System.out.println("Reading MSA"); Sequences alignment =
	 * FastaRW.readFastaAlign(msaFile); SubstitutionMatrix Q =
	 * SMD.create(substitutionModel); MSAData msa = new
	 * MSAData(Q.getSequenceType(),alignment);
	 * 
	 * 
	 * 
	 * System.out.println("Constructing initial guest tree"); double [][] distances
	 * = DistanceMatrix.compute(msa,Q, DistanceMatrix.KIMURA); Pair<Tree , double []
	 * > gTree = NJConstructor.RasmussenNJ(distances,
	 * msa.getSequencesIdentifiers()); gTree.first.show();
	 * 
	 * 
	 * System.out.println(); System.out.println("Sequence Evolution Model");
	 * SubstitutionModel smodel = new SubstitutionModel(msa, Q, gTree, false);
	 * System.out.print("Tree Likelihood: "); smodel.printLogLikelihood();
	 * 
	 * 
	 * 
	 * Pair<Node,Node> nniPair = ProposeTree.proposeRandomNni(gTree.first);
	 * System.out.println("Node_a = " + nniPair.first.id + "\t" + "Node_b = " +
	 * nniPair.second.id ); ProposeTree.performNni(gTree.first, nniPair);
	 * gTree.first.show();
	 * 
	 * System.out.println(); System.out.println("Sequence Evolution Model"); smodel
	 * = new SubstitutionModel(msa, Q, gTree, false);
	 * System.out.print("Tree Likelihood: "); smodel.printLogLikelihood();
	 * 
	 * 
	 * System.out.println("\n"); int [] gene2species =
	 * LeavesMapping.getLeavesMap(gTree.first,sTree,mappingFile);
	 * LinkedHashMap<Integer, Integer> recon = Reconciliation.reconcile(gTree.first,
	 * sTree, gene2species); Reconciliation.labelEvents(gTree.first,recon);
	 * 
	 * long startTime = System.nanoTime();
	 * TopologyPriors.birthDeathTreePriorFull(gTree.first, sTree, recon,params);
	 * long stopTime = System.nanoTime(); System.out.println( (stopTime - startTime)
	 * / Math.pow(10,9));
	 */
	// Testing Area
	// System.out.println("---------------------");
	// System.out.println("Testing Area:");
	// -3457.9427221080273
	// testSystem.validateSeqEvoModelPrime(smodel);
	// testSystem.validateMLE_d(msa,Q);
	// F( R ≤ r ) = 0.4448460821636304
	// testSystem.TestRateModel(params) ;
	// long seed = 1000;
	// Random generator = new Random(seed);
	// double num = generator.nextDouble() * (0.5);
	// System.out.println(num);

	// gtree.show();
	// After adding implied species node show the reconciliation
	/*
	 * String outputPath = "/Users/auwnm/Desktop/Running/spimap-1.1/exe/bin/";
	 * Reconciliation.show(gtree,htree, recon); visTree.showRecon(gtree, htree,
	 * outputPath);
	 */

	// Topology Prior
	// Log form of expression of SubTreeProb for (S>1), (S==1) cases ...
	// subTreeProb = Math.log(phist) + Math.log(p1) + (s-1) * ( Math.log(birthRate)
	// - Math.log(deathRate) + Math.log(p0) ) - (s+1) * Math.log(a);
	// subTreeProb = Math.log(phist) + Math.log(p1) + 2 * Math.log(a);

	/*
	 * LogDouble geneDLModelLikelihood =
	 * TopologyPriors.birthDeathTreePriorFull(gsMap, pi.geneBirthRate,
	 * pi.geneDeathRate, speciesExtinctionProb,log);
	 * pi.geneTree.setArcTimesFromVertexTimes(); log.write("\nSampled time");
	 * pi.geneTree.logTree(null, log);
	 */

	// From main file ...
	// Old parameters setting
	// private static double birthRate = 0.1;
	// private static double deathRate = 0.01;
	// private static double edgeRatePDMean = 0.5;
	// private static double edgeRatePDCV = 1.0;

	// private static String substitutionModel = "JTT";
	// private static String sTreeFile =
	// "/Users/auwnm/Desktop/Running/spimap-1.1/exe/bin/mammals1.stree";
	// //fungi.stree"; //small.stree";
	// private static String msaFile =
	// "/Users/auwnm/Desktop/Running/spimap-1.1/exe/bin/mammals.msa";
	// //100.nt.align"; //small.msa";
	// private static String mappingFile =
	// "/Users/auwnm/Desktop/Running/spimap-1.1/exe/bin/mammals.map";
	// //fungiMap.txt"; //small.map";

	/*
	 * // // // MSA based approach to estimate the rates .... // // // // Rate
	 * Parameters Estimation // Maximum likelihood estimation available at //
	 * http://en.wikipedia.org/wiki/Gamma_distribution public static void
	 * estimateEdgeRates(Parameters pi,Reconciliation geneRecon,Reconciliation []
	 * domainRecon,Log log) throws IOException {
	 * 
	 * StringBuilder logStr = new StringBuilder();
	 * 
	 * // Total arc time of host (species) tree (including stem arc) double totime =
	 * pi.speciesTree.getTotalTime(); assert(totime != 0.0);
	 * 
	 * double rateAtSite; ArrayList<Double> rates; for(int r=0 ; r <
	 * pi.numberOfDomains ; r++) {
	 * 
	 * rates = new ArrayList<Double>(); // Map where patterns (unique columns) are
	 * keys and [first position, count] of the patterns are values.
	 * LinkedHashMap<String, int[]> patterns = pi.domMSA[r].getPatterns();
	 * for(String col: patterns.keySet() ) { int noSubs = 0; int count =
	 * patterns.get(col)[1]; for(int i=0 ; i < col.length()-1 ; i++) if(
	 * col.charAt(i) != col.charAt(i+1) ) noSubs++;
	 * 
	 * rateAtSite = (double) (count*noSubs) / totime ;
	 * 
	 * if(rateAtSite>0.0) rates.add(rateAtSite); }
	 * 
	 * // Step 2: Estimating new rates based on these times double x; double sum =
	 * 0.0; double logSum = 0.0; int N= rates.size(); for(int idx=0; idx < N ; idx++
	 * ) { x = rates.get(idx); sum += x; logSum += Math.log(x); }
	 * 
	 * double s = Math.log( sum ) - Math.log( (double) N ) - ( 1.0 / (double) N ) *
	 * logSum; double k = ( 3.0 - s + Math.sqrt( (s-3.0) * (s-3.0) + 24*s ) ) /
	 * (12.0 * s) ; double theta = ( 1.0 / (k*N) ) * sum;
	 * 
	 * pi.domEdgeRates[r][0] = k*theta; pi.domEdgeRates[r][1] = Math.sqrt(k)/k;
	 * 
	 * logStr.append(pi.domainTree[r].name + " : " + "estimated rate parameters = ["
	 * + pi.domEdgeRates[r][0] + "," + pi.domEdgeRates[r][0] + "]");
	 * logStr.append("\n"); }
	 * 
	 * if(log != null && logStr.length() != 0 ) log.write(logStr.toString()); }
	 * 
	 * 
	 * // // // // Newer ......... // // Rate Parameters Estimation // Maximum
	 * likelihood estimation available at //
	 * http://en.wikipedia.org/wiki/Gamma_distribution public static void
	 * estimateEdgeRates1(Parameters pi,Reconciliation geneRecon,Reconciliation []
	 * domainRecon,Log log) throws IOException {
	 * 
	 * int maxIterations = 5; int totalIterations = 100000; StringBuilder logStr =
	 * new StringBuilder();
	 * 
	 * LinkedHashMap<Integer, double []> domArcTimes = new LinkedHashMap<Integer,
	 * double []>();
	 * 
	 * LogDouble geneDLModelLikelihood,domainDLModelLikelihood,rateModelLikelihoods;
	 * double [] speciesExtinctionProb =
	 * BirthDeathProbs.calcExtinctionProb(pi.speciesTree,pi.geneBDRates); double []
	 * extinctionProb;
	 * 
	 * // Total arc time of host (species) tree (including stem arc) double totime =
	 * pi.speciesTree.getTotalTime(); assert(totime != 0.0);
	 * 
	 * // Initialize rate parameters with random values1 for(int r=0 ; r <
	 * pi.numberOfDomains ; r++) {
	 * 
	 * double [] r8 = new double[pi.domainTree[r].nnodes]; for(int idx=0; idx <
	 * pi.domainTree[r].nnodes ; idx++ ) { r8[idx] = pi.domainTree[r].bl[idx] /
	 * totime ; } Pair<Double,Double> rates = getMLERates(r8); pi.domEdgeRates[r][0]
	 * = rates.first; pi.domEdgeRates[r][1] = rates.second;
	 * 
	 * }
	 * 
	 * 
	 * 
	 * LogDouble initialLikelihood,currentLikelihood,maxLikelihood;
	 * initialLikelihood = new LogDouble(0.0); currentLikelihood = new
	 * LogDouble(0.0);
	 * 
	 * int iterations = 0; while(iterations < maxIterations ) {
	 * 
	 * // Step 1: Pick time sample which has maximum likelihood maxLikelihood = new
	 * LogDouble(0.0); for (int n=0 ; n<totalIterations ; n++ ) { LogDouble
	 * newLikelihood = new LogDouble(1.0); geneDLModelLikelihood =
	 * TopologyPriors.getBirthDeathTreePriorFull(geneRecon,pi.geneBDRates,
	 * speciesExtinctionProb,log); newLikelihood.mult(geneDLModelLikelihood);
	 * for(int r=0 ; r < pi.numberOfDomains ; r++) { extinctionProb =
	 * BirthDeathProbs.calcExtinctionProb(pi.geneTree,pi.domBDRates[r]);
	 * domainDLModelLikelihood =
	 * TopologyPriors.getBirthDeathTreePriorFull(domainRecon[r], pi.domBDRates[r],
	 * extinctionProb,log); rateModelLikelihoods =
	 * RateModel.getLogLikelihood(pi.domEdgeRates[r],pi.domainTree[r],log);
	 * 
	 * newLikelihood.mult(domainDLModelLikelihood);
	 * newLikelihood.mult(rateModelLikelihoods); }
	 * 
	 * 
	 * if(newLikelihood.greaterThan(maxLikelihood)) { for(int r=0 ; r <
	 * pi.numberOfDomains ; r++) { double [] arcTimes =
	 * Arrays.copyOf(pi.domainTree[r].at, pi.domainTree[r].at.length);
	 * domArcTimes.put(r, arcTimes); } maxLikelihood = newLikelihood; }
	 * 
	 * } currentLikelihood = maxLikelihood; if(
	 * currentLikelihood.lessThan(initialLikelihood) ) break;
	 * 
	 * // Step 2: Estimating new rates based on these times for(int r=0 ; r <
	 * pi.numberOfDomains ; r++) { double [] r8 = new
	 * double[pi.domainTree[r].nnodes]; for(int idx=0; idx < pi.domainTree[r].nnodes
	 * ; idx++ ) { if( pi.domainTree[r].bl[idx] <= 0 && domArcTimes.get(r)[idx] ==
	 * 0.0 ) { logStr.
	 * append("\nError: Either Branch length less than or equal to zero or arc time is zero"
	 * ); r8[idx] = 1.0 / pi.domainTree[r].getTotalTime(); } else r8[idx] =
	 * pi.domainTree[r].bl[idx] / domArcTimes.get(r)[idx] ; }
	 * 
	 * Pair<Double,Double> rates = getMLERates(r8); pi.domEdgeRates[r][0] =
	 * rates.first; pi.domEdgeRates[r][1] = rates.second; }
	 * 
	 * logStr.append("\nIteration = " + iterations + "\tInitial likelihood = " +
	 * initialLikelihood.getLogValue() + "\tCurrent Likelihood = " +
	 * currentLikelihood.getLogValue() );
	 * 
	 * 
	 * initialLikelihood = currentLikelihood; iterations++; }
	 * 
	 * 
	 * if(log != null ) log.write(logStr.toString()); else
	 * System.out.println(logStr.toString());
	 * 
	 * }
	 * 
	 * // // Older ......... // // Rate Parameters Estimation // Maximum likelihood
	 * estimation available at // http://en.wikipedia.org/wiki/Gamma_distribution
	 * public static void estimateEdgeRates(Parameters pi,Reconciliation
	 * geneRecon,Reconciliation [] domainRecon,Log log) throws IOException {
	 * 
	 * StringBuilder logStr = new StringBuilder(); double precision = 1e-3; int
	 * totalIterations = 100000; int maxAllowedIterations = 1;
	 * 
	 * int nnodesDom; double [] domArcTimes; double [] extinctionProb; double
	 * oldMean, newMean; LogDouble newLikelihood,oldLikelihood,maxLikelihood;
	 * 
	 * 
	 * // Total arc time of host (species) tree (including stem arc) double totime =
	 * pi.speciesTree.getTotalTime(); assert(totime != 0.0);
	 * 
	 * 
	 * SequenceType st = pi.Q.getSequenceType();
	 * 
	 * double rateAtSite; ArrayList<Double> rates;
	 * 
	 * for(int r=0 ; r < pi.numberOfDomains ; r++) {
	 * 
	 * // Initialize rate parameters with random values pi.domEdgeRatePDMean[r] =
	 * Math.random(); pi.domEdgeRatePDCV [r] = Math.random();
	 * 
	 * int allowedIterations = 0; oldLikelihood = new LogDouble(0.0); maxLikelihood
	 * = new LogDouble(0.0); nnodesDom = pi.domainTree[r].nnodes; domArcTimes = new
	 * double[nnodesDom];
	 * 
	 * oldMean = 0.0; newMean = pi.domEdgeRatePDMean[r];
	 * 
	 * // Please converge while( Math.abs(newMean - oldMean) > precision &&
	 * allowedIterations < maxAllowedIterations ) {
	 * 
	 * // Step 1: Estimate times which has maximum likelihood for (int n=0 ;
	 * n<totalIterations ; n++ ) {
	 * TopologyPriors.birthDeathTreePriorFull(geneRecon,pi.geneBirthRate,pi.
	 * geneDeathRate,speciesExtinctionProb,null); extinctionProb =
	 * BirthDeathProbs.calcExtinctionProb(pi.geneTree,pi.domBirthRate[r],pi.
	 * domDeathRate[r]); TopologyPriors.birthDeathTreePriorFull(domainRecon[r],
	 * pi.domBirthRate[r], pi.domDeathRate[r], extinctionProb,null); newLikelihood =
	 * RateModel.getLogLikelihood(pi.domEdgeRatePDMean[r],pi.domEdgeRatePDCV[r],pi.
	 * domainTree[r],null);
	 * 
	 * if(newLikelihood.greaterThanOrEquals(oldLikelihood)) for(int idx=0 ;
	 * idx<nnodesDom ; idx++) { domArcTimes[idx] =
	 * pi.domainTree[r].nodes.get(idx).at; maxLikelihood = newLikelihood; }
	 * 
	 * oldLikelihood = newLikelihood; } oldMean = pi.domEdgeRatePDMean[r] ;
	 * 
	 * // Step 2: Estimating new rates based on these times int N=0; double x;
	 * double sum = 0.0; double logSum = 0.0; double [] branchLengths =
	 * pi.domLengths.get(r); for(int idx=0; idx < branchLengths.length ; idx++ ) {
	 * 
	 * if( branchLengths[idx] <= 0 && domArcTimes[idx] == 0.0 ) { logStr.
	 * append("Error: Either Branch length less than or equal to zero or arc time is zero"
	 * ); continue; }
	 * 
	 * x = branchLengths[idx] / domArcTimes[idx] ; sum += x; logSum += Math.log(x);
	 * N++; } double s = Math.log( sum ) - Math.log( (double) N ) - ( 1.0 / (double)
	 * N ) * logSum; double k = ( 3.0 - s + Math.sqrt( (s-3.0) * (s-3.0) + 24*s ) )
	 * / (12.0 * s) ; double theta = ( 1.0 / (k*N) ) * sum;
	 * 
	 * pi.domEdgeRatePDMean[r] = k*theta; pi.domEdgeRatePDCV[r] = Math.sqrt(k)/k;
	 * newMean = pi.domEdgeRatePDMean[r];
	 * 
	 * //Completing its iteration allowedIterations++;
	 * 
	 * if(allowedIterations == maxAllowedIterations ) { newLikelihood =
	 * RateModel.getLogLikelihood(pi.domEdgeRatePDMean[r],pi.domEdgeRatePDCV[r],pi.
	 * domainTree[r], pi.domLengths.get(r),null); logStr.append("\n");
	 * logStr.append("MLE Mean = " + newMean + " Likelihood = " +
	 * newLikelihood.getLogValue()); } else { logStr.append("Iteration = " +
	 * allowedIterations + " Mean = " + oldMean + " Likelihood = " +
	 * maxLikelihood.getLogValue()); }
	 * 
	 * } logStr.append("\n"); } if(log != null ) log.write(logStr.toString()); }
	 * 
	 */

	/*
	 * public RealInterval getTimeInterval(Node guestNode) {
	 * 
	 * Integer hostID = this.nodeMap.get(guestNode.id); double time =
	 * this.host.nodes.get(hostID).at;
	 * 
	 * if( !(time>=0.0) ) { System.out.println("Error: Slice time of host node " +
	 * hostID + " is less than or equal to zero."); return null; }
	 * 
	 * RealInterval timeInterval = new RealInterval(0, time, true, true); return
	 * timeInterval; }
	 */

	public static void TestRateModel(Parameters params) {
		System.out.println();
		System.out.println("Rate Model");
		double edgeRatePDMean = 0.5;
		double edgeRatePDCV = 1.0;
		GammaDistribution edgeRatePD = new GammaDistribution(edgeRatePDMean, edgeRatePDCV);
		System.out.println("Gamma Distribution Mean = " + edgeRatePD.getMean() + " Cv = " + edgeRatePD.getCV());

		double length = 0.1;
		double time = 0.33984136668700426;
		// F( R ≤ r ) = 0.4448460821636304

		double prob = edgeRatePD.getCDF(length / time);
		System.out.println("l = " + length + " t = " + time + " r = " + length / time + " F( R ≤ r ) = " + prob);
	}

	public static void checkEigenDecomposition() {
		System.out.println();

		DenseMatrix64F A = new DenseMatrix64F(2, 2);
		DenseMatrix64F rA = new DenseMatrix64F(2, 2);

		DenseMatrix64F E = new DenseMatrix64F(4, 1);
		DenseMatrix64F V = new DenseMatrix64F(2, 2);
		DenseMatrix64F iV = new DenseMatrix64F(2, 2);
		// DenseMatrix64F td = new DenseMatrix64F(4, 1);
		DenseMatrix64F temp = new DenseMatrix64F(2, 2);

		double[] data = { 2, 3, 2, 1 };
		A.setData(data);

		// System.out.println(a.get(1, 1));

		System.out.println("A Matrix =");
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				System.out.print(String.format("%.2f", A.get(i, j)) + "\t");
			}
			System.out.println();
		}

		EigenDecomposition<DenseMatrix64F> eigFact = DecompositionFactory.eigGeneral(2, true);
		if (!eigFact.decompose(A)) {
			throw new RuntimeException("Unable to decompose eigensystem for substitution model.");
		}
		AdditionalEJMLOps.getEigensystemSolution(2, eigFact, E, V);
		CommonOps.invert(V, iV);

		// AdditionalEJMLOps.elementExp(2, E, 1, td);
		AdditionalEJMLOps.multDiagA(2, E, iV, temp);
		CommonOps.mult(V, temp, rA);

		System.out.println("A Matrix =");
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				System.out.print(String.format("%.2f", rA.get(i, j)) + "\t");
			}
			System.out.println();
		}

	}

	public static void validateSeqEvoModelPrime(SubstitutionModel smodel) throws IOException {
		System.out.println();
		System.out.println("Prime @ Initial Tree Log Likelihood");
		String primeInitTreeFile = "/Users/auwnm/Desktop/Running/spimap-1.1/exe/bin/mammalsPrime.gtree";
		Tree gTree = NewickReader.readNewickTreeFile(primeInitTreeFile);
		smodel.updateTreeParameter(gTree);
		smodel.printLogLikelihood();
	}

	public static void validateMLE_d(MSAData msa, SubstitutionMatrix msm) throws RuntimeException {

		int seq1 = 2;
		int seq2 = 4;

		double N[][] = DistanceMatrix.countReplacements(msa.dataAsStrings[seq1], msa.dataAsStrings[seq2],
				msa.getSequenceType());

		System.out.println("Just to observe the function response ....");
		double delta = 0.005;
		for (double i = 1; (i * delta) < 10.0; i++) {
			System.out.println((i * delta) + "\t" + msm.getSeqDerivLogLikelihood(i * delta, N) + "\t"
					+ msm.getSeqLogLikelihood(i * delta, N));
		}

		System.out.println();
		LikelihoodFunction lFunction = new LikelihoodFunction(msm, N);
		final int maxIterations = 1000;
		final double relativeAccuracy = 1.0e-12;
		final double absoluteAccuracy = 1.0e-8;
		UnivariateSolver nonBracketing = new BrentSolver(relativeAccuracy, absoluteAccuracy);
		double dist = nonBracketing.solve(maxIterations, lFunction, 0.0001, SubstitutionMatrix.MAX_MARKOV_TIME);
		System.out.println("-----------------------------------------");
		System.out.println("Distance approximation  = " + dist);
		System.out.println("Log likelihood Function value @ solution  = " + msm.getSeqLogLikelihood(dist, N));

	}

}
