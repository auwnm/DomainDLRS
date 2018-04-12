package posteriorAnalysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;

import phylogeny.Mapping;
import phylogeny.NewickReader;
import phylogeny.NewickWriter;
import phylogeny.Node;
import phylogeny.Reconciliation;
import phylogeny.Tree;

public class PosteriorMain {

	/*
	 * Map Tree: The use of the term MAP tree has been often misused in the
	 * literature. It has sometimes been used to describe the tree associated with
	 * the sampled state in the MCMC chain that has the highest posterior
	 * probability density. This is problematic, because the sampled state with the
	 * highest posterior probability density may just happen to have extremely good
	 * branch lengths on an otherwise fairly average tree topology. A better
	 * definition of the MAP tree topology is the tree topology that has the
	 * greatest posterior probability, averaged over all branch lengths and
	 * substitution parameter values. For data sets that are very well resolved, or
	 * have a small number of taxa, this is easily calculated by just determining
	 * which tree topology has been sampled the most often in the chain. However for
	 * large data sets it is quite possible that every sampled tree has a unique
	 * topology. In this case, the alternative definition of maximum clade
	 * credibility tree above can be used.
	 * 
	 * Maximum clade credibility tree: A natural candidate for a point estimate is
	 * the tree with the maximum product of the posterior clade probabilities. To
	 * the extent that the posterior probabilities of different clades are additive,
	 * This definition is an estimate of the total probability of the given tree
	 * topology i.e. it provides a way of estimating the maximum a posteriori tree
	 * topology.
	 * 
	 * A maximum clade credibility tree is a tree that summarises the results of a
	 * Bayesian phylogenetic inference. Whereas a majority-rule tree combines the
	 * most common clades, and usually yields a tree that wasn't sampled in the
	 * analysis. The maximum-credibility method evaluates each of the sampled
	 * posterior trees. Each clade within the tree is given a score based on the
	 * fraction of times that it appears in the set of sampled posterior trees, and
	 * the product of these scores are taken as the tree's score. The tree with the
	 * highest score is then the maximum clade credibility tree.
	 * 
	 * 95% credible set of tree topologies: For data sets that are very well
	 * resolved, you might consider reporting not a single tree, but the smallest
	 * set of all tree topologies that accounts for 95% of the posterior
	 * probability. This is known as the 95% credible set of tree topologies.
	 * 
	 * References: Important points taken from the BEAST
	 * http://beast.bio.ed.ac.uk/summarizing-posterior-trees
	 * https://en.wikipedia.org/wiki/Maximum_clade_credibility_tree (Drummond,
	 * Alexei; Rambaut, Andrew. "Summarizing posterior trees". Wikipedia)
	 */

	// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesRun1/";
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results4/ZNF468/";
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/mccTrees/mccZNF468Tree.txt";
	// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesMCC2/mccZNFPRDM9Tree.txt"
	// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesGeneRun/";
	// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesGeneMap/mapGenePRDM9Tree.txt";
	// "/Users/auwnm/Desktop/Running/Results/jPrime/jPrimeGeneRun/"
	// "/Users/auwnm/Desktop/Running/Results/jPrime/jPrimeGeneMap/mapGenePRDM9.txt"

	public static String familyPath;
	public static String mccFileName;

	public static String postAnalysisName;
	public static String binTreeFileName = familyPath + "bintree";
	public static String simOutFileName;

	// Reconciliation variables
	public static boolean searchAlternatives = false;
	public static Tree hostTree;
	public static Mapping guest2host;
	public static Reconciliation recon;

	public static void main(String[] args) throws IOException {

		// computeRFDistances1() ;
		// readDomainDLRSSimRuns();
		// readMrBayesSimRuns();

		boolean computeStatistics = false;
		boolean domainDLRSPostAnalysis = true;
		boolean mrBayesPostAnalysis = false;
		boolean jPrimePostAnalysis = false;

		if (computeStatistics)
			// computeStatis();
			computeStatis1();

		if (domainDLRSPostAnalysis)
			PerformDomainDLRSPostAnalysis();

		if (mrBayesPostAnalysis)
			PerformMrBayesPostAnalysis();

		if (jPrimePostAnalysis)
			PerformJPrimePostAnalysis();

	}

	public static void PerformDomainDLRSPostAnalysis() throws IOException {

		postAnalysisName = familyPath + "postAnalysis.txt";

		// familyPath + "ZNF91_dom0.mcmc";
		// familyPath + "ZNF91_dom1.mcmc";
		// familyPath + "ZNF91_gene.mcmc";

		familyPath = "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results6/Run_PRDM9/FamPRDM9_Primates/"; // PRDM9_Mouse_Small_LargeRun/";
																												// //PRDM9_Manual_Removed_Identical/";
																												// //PRDM9_Mouse_Small/";
																												// //ZNF91_Manual/";
																												// //ZNF91/";
		mccFileName = "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results6/Run_PRDM9/FamPRDM9_Primates/mccPRDM9_Primates_try1.txt"; // mccZNF91Tree_alt_try.txt";
																																			// //
																																			// AlternativesMapZNF/mccPRDM9ZNF_Manual_Removed_Identical.txt";
																																			// //mccZNF91Tree_Manual.txt";
																																			// //mapZNF91Tree1.txt";

		// familyPath =
		// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results7/ZNF468_topology_prt/";
		// mccFileName =
		// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results7/ZNF468_topology_prt/mccZNF468Tree.txt";

		String runFile = familyPath + "PRDM9_dom2.mcmc"; // "ZNF468_dom1.mcmc"; //"ZNF91_dom1.mcmc";
															// //"PRDM9_dom2.mcmc"; //"ZNF91_dom1.mcmc";

		boolean header = true;
		boolean domainPostAnalysis = true;
		ArrayList<Sample> samples;
		if (domainPostAnalysis)
			samples = PosteriorReader.readDomainRunFile(runFile, header);
		else
			samples = PosteriorReader.readGeneRunFile(runFile, header);

		double domainDLRSBurnin = 0.25;
		double domainDLRSCredibleSet = 0.95;

		// Reading reconciliation settings
		if (searchAlternatives) {
			String hostTreeFile = familyPath + "mapGeneTree.txt";
			String mappingFile = familyPath + "ZNF.map";

			guest2host = new Mapping();
			guest2host.readMapping(mappingFile);
			hostTree = NewickReader.readNewickTreeFile(hostTreeFile);
			hostTree.setName("hostTree");
			hostTree.setDepths(hostTree.root, 0);
			hostTree.setTimesFromOriginalLengths();
		}

		runAnalysis(samples, domainDLRSBurnin, domainDLRSCredibleSet);
	}

	public static void PerformMrBayesPostAnalysis() throws IOException {

		/*
		 * //For Mr Bayes posterior //Reading Single Column file String
		 * singleColFileName =
		 * "/Users/auwnm/Desktop/Running/Results/MrBayesResults/PRDM9.mbpost.out";
		 * postAnalysisName = singleColFileName + ".nwk"; boolean singleColHeader =
		 * false; double singleBurnin = 0.0; ArrayList<Sample> singleColSamples =
		 * readSingleColFile(singleColFileName,singleColHeader);
		 * runAnalysis(singleColSamples, singleBurnin, 1.0);
		 */

		familyPath = "/Users/auwnm/Desktop/Running/Results/MrBayes/PRDM9ZNFRun_RI/"; // PRDM9ZNFRun/";
																						// //MRBayesGeneRun/";
		mccFileName = "/Users/auwnm/Desktop/Running/Results/MrBayes/PRDM9ZNFRun_RI/mccPRDM9ZNFTree_Primates_RI.txt"; // PRDM9ZNFRun/mapPRDM9.mcc";
																														// //MRBayesGeneMap/mapGenePRDM9Tree.txt";

		boolean onlyMapTree = false;
		ArrayList<String> treeStrList = PosteriorReader
				.readMrBayesPosterior(familyPath + "PRDM9_Primates_ZNF_RI.nex.trprobs", onlyMapTree); // "PRDM9_Mouse_Small_ZNF_RI.nex.trprobs"
																										// //
																										// "PRDM9_Mouse_Small_ZNF.nex.trprobs"
																										// "PRDM9_Mouse_Gene.nex.trprobs"
																										// //"PRDM9_NonIdentical.nex.trprobs"
																										// //"PRDM9.nex.trprobs"
		ArrayList<Sample> samples = new ArrayList<Sample>();
		if (!onlyMapTree) {
			for (int j = 0; j < treeStrList.size(); j++) {
				Sample sample = new Sample();
				sample.iteration = j;
				sample.tree = treeStrList.get(j).trim();
				samples.add(sample);
			}
			double mrBayesBurnin = 0.0;
			double mrBayesCredibleSet = 1.0;
			runAnalysis(samples, mrBayesBurnin, mrBayesCredibleSet);
		} else {

			String mapTreeStr = treeStrList.get(0);
			File outputFile = new File(mccFileName);
			BufferedWriter bwOut = new BufferedWriter(new FileWriter(outputFile, true));
			bwOut.write(mapTreeStr);
			bwOut.newLine();
			bwOut.close();

		}

	}

	public static void PerformJPrimePostAnalysis() throws IOException {

		familyPath = "/Users/auwnm/Desktop/Running/Results/jPrime/jPrimeGeneRun/";
		mccFileName = "/Users/auwnm/Desktop/Running/Results/jPrime/jPrimeGeneMap/mapGene611.txt";

		boolean onlyMapTree = false;
		ArrayList<String> treeStrList = PosteriorReader.readJPrimePosterior(familyPath + "ZNF611.mcmc", onlyMapTree);
		ArrayList<Sample> samples = new ArrayList<Sample>();
		for (int j = 0; j < treeStrList.size(); j++) {
			Sample sample = new Sample();
			sample.iteration = j;
			sample.tree = treeStrList.get(j).trim();
			samples.add(sample);
		}
		double jPrimeBurnin = 0.25;
		double jPrimeCredibleSet = 0.95;
		runAnalysis(samples, jPrimeBurnin, jPrimeCredibleSet);

	}

	// DomainDLRS
	// ZNF91 case:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results6/Run_2m/Fam91_paper/mccZNF91Tree_alt.txt";
	// ZNF558 case:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/mccZNFTrees/mccZNF558Tree.txt";
	// ZNF468 case:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/mccZNFTrees/mccZNF468Tree.txt";
	// ZNF611 case:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/mccZNFTrees/mccZNF611Tree_new.txt";
	// // Old is actually new
	// ZNF679 case:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/mccZNFTrees/mccZNF679Tree.txt";
	// ZNF764 case:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/mccZNFTrees/mccZNF764Tree.txt";
	// PRDM9 Primates:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results6/Run_PRDM9/FamPRDM9_Primates/mccPRDM9_Primates.txt";
	// PRDM9 Mouse:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results5/PRDM9_Mouse_Small_LargeRun/mccPRDM9_Mouse_Small_LargeRun.txt";

	public static void computeStatis() throws IOException {

		// domainDLRS
		// String mccFileName =
		// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results6/Run_PRDM9/FamPRDM9_Primates/mccPRDM9_Primates.txt";
		// //"/Users/auwnm/Desktop/Running/Results/DomainDLRS/results5/PRDM9_Mouse_Small_LargeRun/mccPRDM9_Mouse_Small_LargeRun.txt";
		// //mccZNFTrees/mccZNF91Tree.txt"; //mccZNF468Tree.txt";
		// String mapGeneTreeFile =
		// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/mapGeneTrees/mapGenePRDM9Tee.txt";
		// //"/Users/auwnm/Desktop/Running/Results/DomainDLRS/mapGeneTrees/mapGenePRDM9Tree_Mouse.txt";
		// //mapGeneTrees/mapGene91Tee.txt";
		// String mappingFile =
		// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results6/Run_PRDM9/FamPRDM9_Primates/ZNF.map";
		// //"/Users/auwnm/Desktop/Running/Results/DomainDLRS/results5/PRDM9_Mouse_Small_LargeRun/ZNF.map";
		// // results4/ZNF91/ZNF.map";

		// mrBayes
		// String mccFileName =
		// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/mapZNF558_MB.txt";//"/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/Mouse_Small_mccPRDM9ZNFTree_rt.txt";
		// String mapGeneTreeFile =
		// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/mapGene558Tee_domainDLRS.txt";//"/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/Mouse_Small_mapGenePRDM9Tree.txt";
		// String mappingFile =
		// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/ZNF558.map";//"/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/Mouse_Small_ZNFPRDM9.map";

		// mrBayes
		String mccFileName = "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/mapZNF611_MB.txt";// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/Mouse_Small_mccPRDM9ZNFTree_rt.txt";
		String mapGeneTreeFile = "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/mapGene611Tee_domainDLRS.txt";// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/Mouse_Small_mapGenePRDM9Tree.txt";
		String mappingFile = "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/ZNF611.map";// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/Mouse_Small_ZNFPRDM9.map";

		// Reading host tree file
		Tree hostTree = NewickReader.readNewickTreeFile(mapGeneTreeFile);
		hostTree.setName("hostTree");
		hostTree.setDepths(hostTree.root, 0);
		hostTree.setTimesFromOriginalLengths();

		// Reading mapping file
		Mapping guest2host = new Mapping();
		guest2host.readMapping(mappingFile);

		// Reading guest tree
		Tree guestTree = NewickReader.readNewickTreeFile(mccFileName);
		guestTree.setName("guestTree");
		guestTree.setDepths(guestTree.root, 0);
		guestTree.setBranchLengthsFromOriginalLengths();

		// Constructing reconciliation
		Reconciliation recon = new Reconciliation(guestTree, hostTree, guest2host);

		int dups = 0;
		int old_dups = 0;
		int recent_dups = 0;

		int bifurcations = 0;
		int losses = recon.noImpliedNodes;
		int totalNodes = guestTree.nnodes;

		double[] supportVector = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
		double[] recentCount = new double[supportVector.length];
		double[] ancientCount = new double[supportVector.length];

		// Collect its recent edges
		for (Node gnode : recon.guest.nodes) {

			double clProb = gnode.originalLen;
			if (recon.nodeMap[gnode.id] != hostTree.root.id) {
				for (int m = 0; m < supportVector.length; m++)
					if (clProb <= supportVector[m])
						recentCount[m]++;
			} else {
				for (int m = 0; m < supportVector.length; m++)
					if (clProb <= supportVector[m])
						ancientCount[m]++;
			}

			if (recon.eventMap[gnode.id] == Reconciliation.EVENT_DUPL) {
				dups++;

				if (recon.nodeMap[gnode.id] == hostTree.root.id)
					old_dups++;
				else
					recent_dups++;
			}

			if (recon.eventMap[gnode.id] == Reconciliation.EVENT_SPEC)
				bifurcations++;

		}

		System.out.println("Number of leaves = " + guestTree.getNumberOfLeaves());

		System.out.println("Number of recent   duplications = " + recent_dups);
		System.out.println("Number of pre root duplications = " + old_dups);
		System.out.println("Number of duplications = " + dups);

		System.out.println("Number of gene induced domain bifurcations = " + bifurcations);
		System.out.println("Number of losses = " + losses);

		System.out.println("Total Number of nodes = " + totalNodes);
		System.out.println("Recent Edges Posterior Distribution:");
		for (int m = 0; m < recentCount.length; m++) {
			System.out.print((int) recentCount[m] + "\t");
		}

		System.out.println("\n");
		System.out.println("Ancient Edges Posterior Distribution:");
		for (int m = 0; m < ancientCount.length; m++) {
			System.out.print((int) ancientCount[m] + "\t");
		}
		System.out.println("\n");

	}

	/////////////////////////////////////////////////////////////////////////////////////////////

	// DomainDLRS
	// ZNF91 case:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results6/Run_2m/Fam91_paper/mccZNF91Tree_alt.txt";
	// ZNF558 case:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/mccZNFTrees/mccZNF558Tree.txt";
	// ZNF468 case:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/mccZNFTrees/mccZNF468Tree.txt";
	// ZNF611 case:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/mccZNFTrees/mccZNF611Tree_new.txt";
	// ZNF679 case:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/mccZNFTrees/mccZNF679Tree.txt";
	// ZNF764 case:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/mccZNFTrees/mccZNF764Tree.txt";
	// PRDM9 Primates:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results6/Run_PRDM9/FamPRDM9_Primates/mccPRDM9_Primates.txt";
	// PRDM9 Mouse:
	// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results5/PRDM9_Mouse_Small_LargeRun/mccPRDM9_Mouse_Small_LargeRun.txt";

	// MrBayes MPR
	// ZNF91 case:
	// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/mapZNF91_MB.txt";
	// // Same gene tree
	// ZNF558 case:
	// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/mapZNF558_MB.txt";
	// // *DomainDLRS mapGeneTreeFile =
	// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/mapGene558Tee_domainDLRS.txt";
	// ZNF468 case:
	// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/mapZNF468_MB.txt";
	// // Same gene tree
	// ZNF611 case:
	// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/mapZNF611_MB.txt";
	// ZNF679 case:
	// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/mapZNF679_MB.txt";
	// ZNF764 case:
	// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/mapZNF764_MB.txt";
	// PRDM9 Primates
	// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/Primate_mccPRDM9_MB.txt"
	// PRDM9 Primates
	// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/Mouse_Small_mccPRDM9ZNFTree_rt.txt";

	public static void computeStatis1() throws IOException {

		// domainDLRS
		// String mccFileName =
		// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/mccZNFTrees/mccZNF611Tree_new.txt";
		// //"/Users/auwnm/Desktop/Running/Results/DomainDLRS/results5/PRDM9_Mouse_Small_LargeRun/mccPRDM9_Mouse_Small_LargeRun.txt";
		// //mccZNFTrees/mccZNF91Tree.txt"; //mccZNF468Tree.txt";
		// String mapGeneTreeFile =
		// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/mapGeneTrees/mapGene611Tee.txt";
		// //"/Users/auwnm/Desktop/Running/Results/DomainDLRS/mapGeneTrees/mapGenePRDM9Tree_Mouse.txt";
		// //mapGeneTrees/mapGene91Tee.txt";
		// String mappingFile =
		// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results4/ZNF611/ZNF.map";
		// //"/Users/auwnm/Desktop/Running/Results/DomainDLRS/results5/PRDM9_Mouse_Small_LargeRun/ZNF.map";
		// // results4/ZNF91/ZNF.map";

		// domainDLRS
		// String mccFileName =
		// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results5/PRDM9_Mouse_Small_LargeRun/mccPRDM9_Mouse_Small_LargeRun.txt";
		// String mapGeneTreeFile =
		// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results5/PRDM9_Mouse_Small_LargeRun/mapGeneTree.txt";
		// String mappingFile =
		// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results5/PRDM9_Mouse_Small_LargeRun/ZNF.map";

		// mrBayes
		// String mccFileName =
		// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/mapZNF679_MB.txt";
		// //Mouse_Small_mccPRDM9ZNFTree_rt.txt";
		// String mapGeneTreeFile =
		// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/mapGene679Tree_rt.txt";
		// //Mouse_Small_mapGenePRDM9Tree.txt";
		// String mappingFile =
		// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/ZNF679.map";
		// //Mouse_Small_ZNFPRDM9.map";

		// mrBayes
		String mccFileName = "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/mapZNF611_MB.txt";// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/Mouse_Small_mccPRDM9ZNFTree_rt.txt";
		String mapGeneTreeFile = "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/mapGene611Tee_domainDLRS.txt";// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/Mouse_Small_mapGenePRDM9Tree.txt";
		String mappingFile = "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/ZNF611.map";// "/Users/auwnm/Desktop/Running/Results/MrBayes/MRBayesZNFMCC2_Rooted/Mouse_Small_ZNFPRDM9.map";

		// Reading host tree file
		Tree hostTree = NewickReader.readNewickTreeFile(mapGeneTreeFile);
		hostTree.setName("hostTree");
		hostTree.setDepths(hostTree.root, 0);
		hostTree.setTimesFromOriginalLengths();

		// Reading mapping file
		Mapping guest2host = new Mapping();
		guest2host.readMapping(mappingFile);

		// Reading guest tree
		Tree guestTree = NewickReader.readNewickTreeFile(mccFileName);
		guestTree.setName("guestTree");
		guestTree.setDepths(guestTree.root, 0);
		guestTree.setBranchLengthsFromOriginalLengths();

		// Constructing reconciliation
		Reconciliation recon = new Reconciliation(guestTree, hostTree, guest2host);

		int dups = 0;
		int bifurcations = 0;
		int losses = 0;

		int[] ndomains = new int[hostTree.nnodes];
		int[] ndups = new int[hostTree.nnodes];
		int[] nlosses = new int[hostTree.nnodes];

		// Collect its recent edges
		for (Node gnode : recon.guestPrime.nodes) {

			Node snode = recon.host.getNode(recon.nodeMapPrime[gnode.id]);

			if (recon.eventMapPrime[gnode.id] == Reconciliation.EVENT_DUPL) {
				ndups[snode.id]++;
				dups++;
			} else {
				ndomains[snode.id]++;
				bifurcations++;
			}

			if (gnode.isImplied) {
				int impliedNodeChildID = gnode.children.get(0).id;
				int impliedID = recon.nodeMapPrime[impliedNodeChildID];

				for (int i = 0; i < snode.nchildren; i++) {

					if (snode.children.get(i).id != impliedID)
						nlosses[snode.children.get(i).id]++;

					losses++;
				}
			}
		}

		assert (losses == recon.noImpliedNodes);

		System.out.println("Total Number of nodes = " + guestTree.nnodes);
		System.out.println("Number of leaves = " + guestTree.getNumberOfLeaves());
		System.out.println("Number of duplications = " + dups);
		System.out.println("Number of gene induced domain bifurcations = " + bifurcations);
		System.out.println("Number of losses = " + losses);

		System.out.println("\n");
		System.out.println("Domain Counts:");

		ArrayList<Node> nodes = hostTree.getTreePostList();

		String[] nodeLabels = getNewickLebels(hostTree);

		// String.format("%-50s",nodeLabels[snode.id])
		for (Node snode : nodes) {
			System.out.print(snode.id + "\t" + nodeLabels[snode.id] + "\t" + ndomains[snode.id] + "\t"
					+ ndups[snode.id] + "\t" + nlosses[snode.id] + "\n");
		}

		System.out.println("\n");

	}

	// write out the newick notation of a tree
	public static void getNewickNodeName(String[] nodeLabels, Node node, int depth) {

		StringBuilder nwkStr = new StringBuilder();
		if (node.nchildren == 0) {
			nwkStr.append(String.format("%s", node.name));
		} else {

			ArrayList<String> childLabels = new ArrayList<String>();
			for (int i = 0; i < node.nchildren; i++) {
				getNewickNodeName(nodeLabels, node.children.get(i), depth + 1);
				childLabels.add(nodeLabels[node.children.get(i).id]);
			}
			Collections.sort(childLabels);
			nwkStr.append("(");

			for (int i = 0; i < childLabels.size() - 1; i++) {
				nwkStr.append(childLabels.get(i));
				nwkStr.append(",");
			}

			nwkStr.append(childLabels.get(childLabels.size() - 1));
			nwkStr.append(")");

		}
		nodeLabels[node.id] = nwkStr.toString();
	}

	// Output the newick notation of a tree in string
	public static String[] getNewickLebels(Tree tree) {
		String[] nodeLabels = new String[tree.nnodes];
		getNewickNodeName(nodeLabels, tree.root, 0);
		return nodeLabels;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////

	public static void computeRFDistances1() throws IOException {
		// String tree1File =
		// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results5/PRDM9_Mouse_Small_LargeRun/mccPRDM9_Mouse_Small_LargeRun.txt";//AlternativesMapZNF/mccPRDM9Tree_Manual.txt";
		// String tree2File =
		// "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results5/PRDM9_Mouse_Small/mccPRDM9ZNFMouseSmall.txt";
		// //AlternativesMapZNF/mccPRDM9K1Tree.txt";

		String tree1File = "/Users/auwnm/Desktop/Running/Results/DomainDLRS/results7/ZNF468_topology_prt/mccZNF468Tree.txt";
		String tree2File = "/Users/auwnm/Desktop/Running/Results/DomainDLRS/mccZNFTrees/mccZNF468Tree.txt";

		Tree tree1 = NewickReader.readNewickTreeFile(tree1File);
		Tree tree2 = NewickReader.readNewickTreeFile(tree2File);

		double rfDist = RobinsonFouldsDistance(tree1, tree2);
		System.out.println(tree1.getNumberOfLeaves() + "\t" + tree2.getNumberOfLeaves() + "\t" + rfDist);

	}

	public static void computeRFDistances() throws IOException {

		// Specify files
		String inferredTrees = "/Users/auwnm/Desktop/Running/Results/Simulations2/ZNF91/simZNF91_MrBayes_ZNF.map";
		String groundTruthPath = "/Users/auwnm/Documents/Jworkspace/JPrIME_runs/GenPhyloData/Simulator/simZNF91/";

		// First reading the inferred trees
		String line;
		BufferedReader br = new BufferedReader(new FileReader(inferredTrees));
		ArrayList<String> inferredTreeStr = new ArrayList<String>();
		while ((line = br.readLine()) != null) {
			if (line.length() == 0)
				break;
			inferredTreeStr.add(line.trim());
		}
		br.close();

		// Writing output to file ...
		File outputFile = new File(inferredTrees + ".rfd");
		BufferedWriter bwOut = new BufferedWriter(new FileWriter(outputFile, true));
		for (int i = 0; i < inferredTreeStr.size(); i++) {

			String groundTruthFile = groundTruthPath + "domainTree" + (i + 1) + "_relaxed.nwk";

			// String groundTruthFile = groundTruthPath + "geneTree_relaxed.nwk";

			Tree groundTruthTree = NewickReader.readNewickTreeFile(groundTruthFile);
			Tree inferedTree = NewickReader.readNewickTreeStr(inferredTreeStr.get(i));

			double rfDist = RobinsonFouldsDistance(groundTruthTree, inferedTree);

			bwOut.write(String.valueOf(rfDist));
			bwOut.newLine();
		}
		bwOut.close();
	}

	public static void readMrBayesSimRuns() throws IOException {

		boolean onlyMapTree;
		String resultPath = "/Users/auwnm/Desktop/Running/Results/Simulations2/ZNF679/simZNF679_MrBayes/";
		simOutFileName = "/Users/auwnm/Desktop/Running/Results/Simulations2/ZNF679/simZNF679_MrBayes_ZNF.map";
		int numberOfFamilies = 100;

		boolean domainRun = true;
		ArrayList<String> treeStrList;

		for (int i = 0; i < numberOfFamilies; i++) {

			if (domainRun) {
				onlyMapTree = true;
				String fileName = resultPath + "domainTree" + (i + 1) + ".nex.trprobs";
				treeStrList = PosteriorReader.readMrBayesPosterior(fileName, onlyMapTree);

				// Jani bemani
				/*
				 * ArrayList<Sample> samples = new ArrayList<Sample>(); for(int j=0 ; j<
				 * treeStrList.size() ; j++) { Sample sample = new Sample(); sample.iteration =
				 * j; sample.tree = treeStrList.get(j).trim(); samples.add(sample); } double
				 * mrBayesBurnin = 0.0; runAnalysis(samples, mrBayesBurnin, 1.0);
				 */

				String mapTreeStr = treeStrList.get(0);

				File outputFile = new File(simOutFileName);
				BufferedWriter bwOut = new BufferedWriter(new FileWriter(outputFile, true));
				bwOut.write(mapTreeStr);
				bwOut.newLine();
				bwOut.close();

			} else {
				onlyMapTree = true;
				String fileName = resultPath + "family" + (i + 1) + ".nex.trprobs";
				treeStrList = PosteriorReader.readMrBayesPosterior(fileName, onlyMapTree);
				String mapTreeStr = treeStrList.get(0);

				File outputFile = new File(simOutFileName);
				BufferedWriter bwOut = new BufferedWriter(new FileWriter(outputFile, true));
				bwOut.write(mapTreeStr);
				bwOut.newLine();
				bwOut.close();
			}

		}

	}

	public static void readDomainDLRSSimRuns() throws IOException {

		boolean HILLCLIMBING = false;
		int numberOfFamilies = 100;
		int numberOfDomains = 2;
		String resultPath = "/Users/auwnm/Desktop/Running/Results/Simulations2/ZNF679/simZNF679_MCMC/";
		simOutFileName = "/Users/auwnm/Desktop/Running/Results/Simulations2/ZNF679/domainDLRS_MCMC_Gene.map";
		int r = 1;

		boolean domainRun = false;
		String runOf;
		if (domainRun)
			runOf = "_dom" + r + ".mcmc";
		else
			runOf = "_gene.mcmc";

		boolean header = true;
		ArrayList<Sample> samples = null;
		for (int i = 0; i < numberOfFamilies; i++) {

			if (HILLCLIMBING) {

				String hillClimbingOutFile = resultPath + "family" + (i + 1) + ".inf";
				String[] data = PosteriorReader.readHCFile(hillClimbingOutFile, numberOfDomains);

				String treeStr;
				if (domainRun)
					treeStr = data[r + 1];
				else
					treeStr = data[0];

				File outputFile = new File(simOutFileName);
				BufferedWriter bwOut = new BufferedWriter(new FileWriter(outputFile, true));
				bwOut.write(treeStr);
				bwOut.newLine();
				bwOut.close();

			} else {

				String mcmcRunFile = resultPath + "family" + (i + 1) + runOf;
				System.out.println("Running analysis on " + mcmcRunFile);

				if (domainRun)
					samples = PosteriorReader.readDomainRunFile(mcmcRunFile, header);
				else
					samples = PosteriorReader.readGeneRunFile(mcmcRunFile, header);

				double domainDLRSBurnin = 0.25;
				double domainDLRSCredibleSet = 0.95;
				runAnalysis(samples, domainDLRSBurnin, domainDLRSCredibleSet);
			}

		}

	}

	// This method will run analysis on output of domainDLRS
	public static void runAnalysis(ArrayList<Sample> samples, double burnin, double sizeOfCredibleSet)
			throws IOException {

		int startIndex = (int) (samples.size() * burnin);

		// Leaf map, split Maps and frequency distribution over tree parameter
		HashMap<String, Integer> name2id = new HashMap<String, Integer>();
		HashMap<Integer, String> id2name = new HashMap<Integer, String>();
		LinkedHashMap<Integer, HashSet<List<Integer>>> splitMaps = new LinkedHashMap<Integer, HashSet<List<Integer>>>();
		LinkedHashMap<Integer, Integer> treeDist = new LinkedHashMap<Integer, Integer>();

		// Constructing trees and converting them into split maps
		for (int i = startIndex; i < samples.size(); i++) {

			String treeStr = samples.get(i).tree;
			Tree tree = NewickReader.readNewickTreeStr(treeStr);

			if (i == startIndex)
				setLeafMaps(tree.getLeafNames(true), name2id, id2name);

			HashSet<List<Integer>> splitsSet = getDescendantsSet(tree, name2id, id2name);
			splitMaps.put(i, splitsSet);
		}

		// Show names to id mapping
		for (Integer key : id2name.keySet())
			System.out.println(key + "\t" + id2name.get(key));

		System.out.println("\n Size of mapping = " + id2name.size());

		// Making frequency distribution over tree parameter
		int totalCases = 0;
		for (int i = startIndex; i < samples.size(); i++) {

			totalCases++;
			HashSet<List<Integer>> splitsSet = splitMaps.get(i);
			if (treeDist.size() == 0) {
				treeDist.put(i, 1);
				continue;
			}

			// Matching with all trees in the distribution
			boolean found = false;
			for (Integer key : treeDist.keySet()) {
				HashSet<List<Integer>> binSet = splitMaps.get(key);
				if (compareDescSet(splitsSet, binSet)) {
					treeDist.put(key, treeDist.get(key) + 1);
					found = true;
					break;
				}
			}
			if (!found)
				treeDist.put(i, 1);

		}

		// Sort the tree distribution and selecting cases say top 95% of the commulative
		// frequency
		LinkedHashMap<Integer, Integer> sortedTreeDist = sortHashMapByValues(treeDist);
		LinkedHashMap<Integer, Integer> truncatedTreeDist = new LinkedHashMap<Integer, Integer>();
		int commulative = 0;
		int limit = (int) (totalCases * sizeOfCredibleSet);
		for (Integer key : sortedTreeDist.keySet()) {

			int frequency = sortedTreeDist.get(key);
			commulative += frequency;

			truncatedTreeDist.put(key, frequency);

			if (commulative > limit)
				break;
		}

		// Show for example top 95% credible set of tree topologies
		boolean showDist = true;
		if (showDist) {
			for (Integer key : truncatedTreeDist.keySet())
				System.out.print(truncatedTreeDist.get(key) + "\n");
			System.out.println("Total cases in tuncated dist. = " + commulative);
		}

		// Writing top bin's trees into files
		boolean writeBinTrees = false;
		int whichBin = 3;
		int binTreeKey = 0;
		if (writeBinTrees) {
			int binKey = 0;
			int binCount = 0;
			for (Integer key : truncatedTreeDist.keySet()) {

				if (binCount >= whichBin)
					break;

				binKey = key;
				binCount++;
			}
			boolean firstTree = true;
			// int caseNo = 0;
			for (int i = startIndex; i < samples.size(); i++) {

				HashSet<List<Integer>> splitsSet = splitMaps.get(i);
				HashSet<List<Integer>> binSet = splitMaps.get(binKey);
				if (compareDescSet(splitsSet, binSet)) {
					// caseNo++;

					if (firstTree) {
						binTreeKey = i;
						firstTree = false;
					}
					/*
					 * //Do not write them, just show them File outputFile = new
					 * File(binTreeFileName + caseNo + ".txt"); BufferedWriter bwOut = new
					 * BufferedWriter(new FileWriter(outputFile)); bwOut.write(samples.get(i).tree);
					 * bwOut.close();
					 */
					// System.out.println(caseNo + ": " + samples.get(i).tree);
				}

			}
		}

		// MCC Part
		// Computing posterior probabilities of different clades
		double totalTruncatedCases = commulative;
		LinkedHashMap<List<Integer>, Double> cladeProbabilities = new LinkedHashMap<List<Integer>, Double>();
		for (Integer key : truncatedTreeDist.keySet()) {

			double fraction = (double) sortedTreeDist.get(key) / totalTruncatedCases;
			HashSet<List<Integer>> splitSet = splitMaps.get(key);

			for (List<Integer> split : splitSet) {
				if (cladeProbabilities.containsKey(split))
					cladeProbabilities.put(split, cladeProbabilities.get(split) + fraction);
				else
					cladeProbabilities.put(split, fraction);
			}
		}

		// Writing clade probabilities
		boolean writeCladeProbabilities = true;
		if (writeCladeProbabilities) {

			int clade_size = 2;
			File outputFile = new File("/Users/auwnm/Desktop/cladeProbs" + clade_size + ".csv");
			BufferedWriter bwOut = new BufferedWriter(new FileWriter(outputFile));

			ArrayList<String> znfCase1 = new ArrayList<String>();
			ArrayList<String> znfCase2 = new ArrayList<String>();
			ArrayList<String> znfCase3 = new ArrayList<String>();
			ArrayList<String> znfCase4 = new ArrayList<String>();
			ArrayList<String> znfCase5 = new ArrayList<String>();

			znfCase1.add("HS_PRDM9_ZNF4");
			znfCase1.add("HS_PRDM9_ZNF7");
			znfCase1.add("HS_PRDM9_ZNF8");

			znfCase2.add("HS_PRDM9_ZNF5");
			znfCase2.add("HS_PRDM9_ZNF6");

			znfCase3.add("HS_PRDM9_ZNF9");
			znfCase3.add("HS_PRDM9_ZNF12");

			znfCase4.add("PT_PRDM9_ZNF5");
			znfCase4.add("PT_PRDM9_ZNF7");
			znfCase4.add("PT_PRDM9_ZNF8");
			znfCase4.add("PT_PRDM9_ZNF13");
			znfCase4.add("PT_PRDM9_ZNF14");

			znfCase5.add("PG_PRDM9_ZNF4");
			znfCase5.add("PG_PRDM9_ZNF6");

			bwOut.write("Clade" + "\t" + "Case" + "\t" + "Species" + "\t" + "Probability" + "\n");
			for (List<Integer> key : cladeProbabilities.keySet()) {

				if (key.size() > clade_size)
					continue;

				boolean spcode_HS = true;
				boolean spcode_PT = true;
				boolean spcode_PG = true;
				boolean spcode_MM = true;

				for (int i = 0; i < key.size(); i++) {
					String speciesCode = id2name.get(key.get(i)).substring(0, 2);

					spcode_HS &= speciesCode.equalsIgnoreCase("HS") ? true : false;
					spcode_PT &= speciesCode.equalsIgnoreCase("PT") ? true : false;
					spcode_PG &= speciesCode.equalsIgnoreCase("PG") ? true : false;
					spcode_MM &= speciesCode.equalsIgnoreCase("MM") ? true : false;

				}

				String speciesStr = "Others";
				if (spcode_HS)
					speciesStr = "Human";
				if (spcode_PT)
					speciesStr = "Chmip";
				if (spcode_PG)
					speciesStr = "Orangutan";
				if (spcode_MM)
					speciesStr = "Macaque";

				int ncases = 0;
				int caseFlag = 0;
				ArrayList<String> znfCase;

				znfCase = znfCase4;
				for (int j = 0; j < znfCase.size(); j++) {

					String znfName = znfCase.get(j);
					Integer id = name2id.get(znfName);

					if (key.contains(id))
						ncases++;

					if (ncases == clade_size)
						break;

				}
				if (ncases == clade_size)
					caseFlag = 4;

				/*
				 * znfCase = znfCase1; for(int j=0 ; j< znfCase.size() ;j++) {
				 * 
				 * String znfName = znfCase.get(j); Integer id = name2id.get( znfName );
				 * 
				 * if(key.contains(id)) ncases ++;
				 * 
				 * if(ncases == clade_size) break;
				 * 
				 * } if(ncases==clade_size) caseFlag = 1;
				 */

				StringBuilder sb = new StringBuilder();
				sb.append("[");

				for (int i = 0; i < key.size(); i++) {
					sb.append(id2name.get(key.get(i)));
					if (i != (key.size() - 1))
						sb.append(",");
				}
				sb.append("]");

				bwOut.write(sb.toString() + "\t" + caseFlag + "\t" + speciesStr + "\t" + cladeProbabilities.get(key));
				bwOut.newLine();
			}
			bwOut.close();
		}

		// Finding MCC tree by scoring each tree based on the product of the posterior
		// clade probabilities
		double product;
		double maxProduct = 0;
		int maxCladeCredibilityTree = 0;
		ArrayList<Integer> alternativeCases = new ArrayList<Integer>();
		int maxAlt = 0;
		double maxProductAlt = 0;
		boolean altCaseExist = true;

		for (Integer key : truncatedTreeDist.keySet()) {
			product = 1.0;
			HashSet<List<Integer>> splitSet = splitMaps.get(key);
			for (List<Integer> split : splitSet)
				product *= cladeProbabilities.get(split);

			if (product > maxProduct) {
				maxCladeCredibilityTree = key;
				maxProduct = product;
			}

			if (searchAlternatives) {
				String guestTreeStr = samples.get(key).tree;
				Tree guestTree = NewickReader.readNewickTreeStr(guestTreeStr);
				recon = new Reconciliation(guestTree, hostTree, guest2host);

				// Collect its recent Duplications
				ArrayList<Node> recentDups = new ArrayList<Node>();
				for (Node gnode : recon.guest.nodes) {
					if (recon.eventMap[gnode.id] == Reconciliation.EVENT_DUPL
							&& recon.nodeMap[gnode.id] != hostTree.root.id) //
						recentDups.add(gnode);
				}

				if (recentDups.size() == 0)
					continue;

				// For ZNF91
				altCaseExist = false;
				boolean condition1 = true;
				int dup_counter_5 = 0;
				for (Node dup : recentDups) {
					for (int n = 0; n < dup.nchildren; n++) {
						Node child = dup.children.get(n);
						if (recon.eventMap[child.id] == Reconciliation.EVENT_DUPL)
							condition1 = false;
					}

					if (recon.nodeMap[dup.id] == 5)
						dup_counter_5++;
				}

				boolean c1 = false, c2 = false, c3 = false;
				for (Node gnode : recon.guest.nodes) {

					if (!gnode.isLeaf() && recon.eventMap[gnode.id] == Reconciliation.EVENT_SPEC) {

						Node child1 = gnode.children.get(0);
						Node child2 = gnode.children.get(1);

						if (child1.isLeaf() && child2.isLeaf()) {
							String str1 = child1.name;
							String str2 = child2.name;

							if (str1.equalsIgnoreCase("HS_G91_ZNF6") && str2.equalsIgnoreCase("PT_G91_ZNF6")
									|| str1.equalsIgnoreCase("PT_G91_ZNF6") && str2.equalsIgnoreCase("HS_G91_ZNF6"))
								c1 = true;

							if (str1.equalsIgnoreCase("HS_G91_ZNF8") && str2.equalsIgnoreCase("PT_G91_ZNF8")
									|| str1.equalsIgnoreCase("PT_G91_ZNF8") && str2.equalsIgnoreCase("HS_G91_ZNF8"))
								c2 = true;

							if (str1.equalsIgnoreCase("HS_G91_ZNF10") && str2.equalsIgnoreCase("PT_G91_ZNF10")
									|| str1.equalsIgnoreCase("PT_G91_ZNF10") && str2.equalsIgnoreCase("HS_G91_ZNF10"))
								c3 = true;
						} else
							continue;
					}

				}

				boolean condition2 = true;
				if (dup_counter_5 > 6)
					condition2 = true;

				if (c1 && c2)
					System.out.println();

				if (condition1 && condition2 && c3)
					altCaseExist = true;

				/*
				 * // For PRDM9 boolean altCaseExist = false;
				 * 
				 * ArrayList<String> list_pt = new ArrayList<String>();
				 * list_pt.add("PT_PRDM9_ZNF5"); list_pt.add("PT_PRDM9_ZNF7");
				 * list_pt.add("PT_PRDM9_ZNF8"); list_pt.add("PT_PRDM9_ZNF13");
				 * list_pt.add("PT_PRDM9_ZNF14");
				 * 
				 * ArrayList<String> list_desc; for(Node dup:recentDups) { list_desc = new
				 * ArrayList<String>(); getDescendants(dup, list_desc);
				 * if(list_pt.containsAll(list_desc) && list_desc.size()==5) altCaseExist =
				 * true; }
				 */

				if (altCaseExist) {
					alternativeCases.add(key);
					if (product > maxProductAlt) {
						maxProductAlt = product;
						maxAlt = key;
					}
				}

			}
		}

		if (altCaseExist)
			System.out.println("Number of alternatives = " + alternativeCases.size());

		boolean showProduct = false;
		if (showProduct) {
			System.out.println("product = " + maxProduct);
			System.out.println("Clade Probabilities:");
			HashSet<List<Integer>> splitSet = splitMaps.get(maxCladeCredibilityTree);
			for (List<Integer> split : splitSet)
				System.out.println(split.toString() + "\t" + String.format("%.2f", cladeProbabilities.get(split)));
		}

		// Labeling maximum clade credibility tree
		if (writeBinTrees)
			maxCladeCredibilityTree = binTreeKey;

		if (searchAlternatives && alternativeCases.size() != 0)
			maxCladeCredibilityTree = maxAlt;// alternativeCases.get(0);

		String mccTreeStr = samples.get(maxCladeCredibilityTree).tree;
		Tree mccTree = NewickReader.readNewickTreeStr(mccTreeStr);
		double[] cladeProbability = new double[mccTree.nnodes];
		for (int i = 0; i < mccTree.nodes.size(); i++) {

			Node node = mccTree.nodes.get(i);
			if (node.nchildren == 0 || node.parent == null) {
				cladeProbability[node.id] = 1.0;
			} else {

				// Get descendants under this node
				ArrayList<String> descendants = new ArrayList<String>();
				getDescendants(node, descendants);

				// Get into the list of integers and sort them
				ArrayList<Integer> list = new ArrayList<Integer>();
				for (String id : descendants)
					list.add(name2id.get(id));
				Collections.sort(list);

				cladeProbability[node.id] = cladeProbabilities.get(list);
			}
		}

		// Show Maximum Clade Credibility tree
		boolean showMCCTree = false;
		if (showMCCTree) {
			System.out.print("MCC Tree:\t");
			String mccStrTemp = NewickWriter.getNewickString(mccTree, cladeProbability, true);
			System.out.println(mccStrTemp);
			System.out.print("MCC Tree:\t");
			System.out.println(mccTreeStr);

			System.out.println();
			for (int l = 0; l < cladeProbability.length; l++)
				System.out.print(cladeProbability[l] + ", ");
		}

		// Write tree(s) into files
		boolean writeMCCTree = true;
		if (writeMCCTree) {
			File outputFile = new File(mccFileName);

			String mccStrTemp = NewickWriter.getNewickString(mccTree, cladeProbability, true);
			BufferedWriter bwOut = new BufferedWriter(new FileWriter(outputFile, true));
			bwOut.write(mccStrTemp);
			bwOut.newLine();
			bwOut.close();
		}

		// Making consensus tree on 95% credible set of trees
		// Finding Most Occuring tree in 95% credible set of trees
		// Finding the best state tree in 95% credible set of trees
		// Writing all analysis to a post analysis file
		boolean makeOtherAnalysis = false;
		if (makeOtherAnalysis) {

			int mostFrequentKey = -1;
			int bestStateKey = -1;
			double maxPosterior = Double.NEGATIVE_INFINITY;
			ArrayList<String> nwkList = new ArrayList<String>();
			for (Integer key : truncatedTreeDist.keySet()) {

				if (mostFrequentKey == -1)
					mostFrequentKey = key;

				if (samples.get(key).overallDensity >= maxPosterior) {
					bestStateKey = key;
					maxPosterior = samples.get(key).overallDensity;
				}

				String nwkStr = samples.get(key).tree;
				nwkList.add(nwkStr);
			}

			String cnsTreeStr = Consensus.getConsensusTree(nwkList);
			String mfrTree = samples.get(mostFrequentKey).tree;
			String bestTree = samples.get(bestStateKey).tree;

			boolean showCNSTree = false; // Show consensus tree
			boolean showMFRTree = false; // Show most frequent tree
			boolean showBSTTree = false; // Show best state Tree

			if (showCNSTree) {
				System.out.print("CNS Tree:\t");
				System.out.println(cnsTreeStr + ";");
			}

			if (showMFRTree) {
				System.out.print("MFR Tree:\t");
				System.out.println(mfrTree);
			}

			if (showBSTTree) {
				System.out.print("BST Tree:\t");
				System.out.println(bestTree);
			}

			StringBuilder postAnalysisStr = new StringBuilder();
			postAnalysisStr.append(">MCC tree with posterior of each clade:" + "\n");
			postAnalysisStr.append(NewickWriter.getNewickString(mccTree, cladeProbability, true) + "\n");
			postAnalysisStr.append(">MCC chain tree with length:" + "\n");
			postAnalysisStr.append(mccTreeStr + "\n");
			postAnalysisStr.append(">Concensus tree" + "\n");
			postAnalysisStr.append(cnsTreeStr + ";" + "\n");
			postAnalysisStr.append(">Most frequent tree" + "\n");
			postAnalysisStr.append(mfrTree + ";" + "\n");
			postAnalysisStr.append(">Best state tree" + "\n");
			postAnalysisStr.append(bestTree + ";" + "\n");

			File outputFile = new File(postAnalysisName);
			BufferedWriter bwOut = new BufferedWriter(new FileWriter(outputFile));
			bwOut.write(postAnalysisStr.toString());
			bwOut.close();
		}

		/*
		 * double maxDupScore = 0; int maxDupScoreTreeKey = -1;
		 * 
		 * double [] supportVector = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}; double
		 * [] supportCount = new double[supportVector.length];
		 * 
		 * double score=0; double lowestPosterior = 1.0; for(Node dup:recentDups) {
		 * 
		 * // Get descendants under this node ArrayList<String> descendants = new
		 * ArrayList<String>(); getDescendants(dup,descendants);
		 * 
		 * // Get into the list of integers and sort them ArrayList<Integer> list = new
		 * ArrayList<Integer>(); for(String id: descendants) list.add(name2id.get(id));
		 * Collections.sort(list);
		 * 
		 * // Get score for this split double clProb = cladeProbabilities.get(list);
		 * 
		 * if(clProb<lowestPosterior) lowestPosterior = clProb;
		 * 
		 * score += clProb; } for(int m=0 ; m<supportVector.length ; m++)
		 * if(lowestPosterior<=supportVector[m]) supportCount[m]++;
		 * System.out.println("Number of recent duplications = " + recentDups.size());
		 * 
		 * 
		 * if(score > maxDupScore) { maxDupScore = score; maxDupScoreTreeKey = key; }
		 * 
		 * for(int m=0 ; m<supportVector.length ; m++){ System.out.print(supportCount[m]
		 * + ","); } System.out.println(" out of posterior bins " +
		 * truncatedTreeDist.size()) ;
		 * 
		 * 
		 */

		/*
		 * // Writing tree dist to file ... boolean writeMCCSet = false; if(writeMCCSet)
		 * {
		 * 
		 * File outputFile = new File(mccOutFileName); BufferedWriter bwOut = new
		 * BufferedWriter(new FileWriter(outputFile)); for(Integer key:
		 * truncatedTreeDist.keySet()) {
		 * 
		 * Tree temp = NewickReader.readNewickTreeStr( samples.get(key).tree ); String
		 * nwStr = NewickWriter.getNewickTopologyString(temp);
		 * 
		 * //Wrtie as many times as frequency suggest (only topology) for(int
		 * i=0;i<truncatedTreeDist.get(key);i++) { bwOut.write(nwStr); bwOut.newLine();
		 * }
		 * 
		 * } bwOut.close();
		 * 
		 * }
		 * 
		 * // Writing samples to file ... File outputFile = new File(mccOutFileName);
		 * BufferedWriter bwOut = new BufferedWriter(new FileWriter(outputFile));
		 * for(int i=startIndex ; i<samples.size() ; i++) {
		 * 
		 * String treeStr = samples.get(i).tree; bwOut.write(treeStr); bwOut.newLine();
		 * 
		 * } bwOut.close();
		 * 
		 */

		/*
		 * // Now just to verify that this maximum clade credibility tree over top 95%
		 * credible set is near to // 1) Consensus tree // 2) Best state tree // 3) Most
		 * Frequent tree
		 * 
		 * boolean doVerify = false; if(doVerify) {
		 * 
		 * // Making consensus on 95% credible set of trees int mostFrequentKey = -1;
		 * int bestStateKey = -1; double maxPosterior = Double.NEGATIVE_INFINITY;
		 * ArrayList<String> nwkList = new ArrayList<String> (); for(Integer key:
		 * truncatedTreeDist.keySet()) {
		 * 
		 * if(mostFrequentKey == -1) mostFrequentKey = key;
		 * 
		 * if(samples.get(key).overallDensity >= maxPosterior) { bestStateKey = key;
		 * maxPosterior = samples.get(key).overallDensity; }
		 * 
		 * 
		 * String nwkStr = samples.get(key).tree; nwkList.add(nwkStr);
		 * 
		 * } System.out.println("Best state tree posterior = " +
		 * samples.get(bestStateKey).overallDensity);
		 * System.out.println("Most frequent tree frequency = " +
		 * truncatedTreeDist.get(mostFrequentKey));
		 * 
		 * String cnsTreeStr = Consensus.getConsensusTree(nwkList); Tree cnsTree =
		 * NewickReader.readNewickTreeStr( cnsTreeStr ); Tree frqTree =
		 * NewickReader.readNewickTreeStr( samples.get(mostFrequentKey).tree ); Tree
		 * bstTree = NewickReader.readNewickTreeStr( samples.get(bestStateKey).tree );
		 * 
		 * 
		 * // Computing Robinson Foulds distance between mcc tree and --- double rfDist1
		 * = RobinsonFouldsDistance(mccTree, cnsTree); double rfDist2 =
		 * RobinsonFouldsDistance(frqTree, cnsTree); double rfDist3 =
		 * RobinsonFouldsDistance(bstTree, cnsTree);
		 * 
		 * System.out.
		 * println("The Robinson Foulds distance between consensus tree and maximum clade credibility tree on 95% credible set of posterior = "
		 * + rfDist1); System.out.
		 * println("The Robinson Foulds distance between consensus tree and most frequent tree on 95% credible set of posterior = "
		 * + rfDist2); System.out.
		 * println("The Robinson Foulds distance between consensus tree and best state tree on 95% credible set of posterior = "
		 * + rfDist3); }
		 */

	}

	// The Robinson–Foulds metric is a way to measure the distance between unrooted
	// phylogenetic trees.
	// It is defined as (A + B) where A is the number of partitions of data implied
	// by the first tree but not the second tree
	// and B is the number of partitions of data implied by the second tree but not
	// the first tree.
	// It is also known as the symmetric difference metric.
	public static double RobinsonFouldsDistance(Tree tree1, Tree tree2) {
		double rfDistance = 0;

		HashMap<String, Integer> name2id = new HashMap<String, Integer>();
		HashMap<Integer, String> id2name = new HashMap<Integer, String>();

		setLeafMaps(tree1.getLeafNames(true), name2id, id2name);

		HashSet<List<Integer>> set1 = getPartitions(tree1, name2id, id2name);
		HashSet<List<Integer>> set2 = getPartitions(tree2, name2id, id2name);

		for (List<Integer> element : set1) {
			if (!set2.contains(element))
				rfDistance += 1.0;
		}
		for (List<Integer> element : set2) {
			if (!set1.contains(element))
				rfDistance += 1.0;
		}

		return rfDistance;
	}

	// Parse tree object into non-trivial partitions and return them in Hashable
	// format
	public static HashSet<List<Integer>> getPartitions(Tree tree, HashMap<String, Integer> name2id,
			HashMap<Integer, String> id2name) {

		HashSet<List<Integer>> splitSet = new HashSet<List<Integer>>();

		for (int i = 0; i < tree.nodes.size(); i++) {

			Node node = tree.nodes.get(i);
			if (node.nchildren == 0 || node.parent == null) {
				continue;
			} else {

				// Get descendants under this node
				ArrayList<String> descendants = new ArrayList<String>();
				getDescendants(node, descendants);

				// Get into the list of integers and sort them
				ArrayList<Integer> list = new ArrayList<Integer>();
				for (String id : descendants)
					list.add(name2id.get(id));
				Collections.sort(list);

				// Making complementary list list U compList = Leaves
				ArrayList<Integer> compList = new ArrayList<Integer>();
				for (Integer id : id2name.keySet()) {
					if (!list.contains(id))
						compList.add(id);
				}
				Collections.sort(compList);

				// Defining split
				ArrayList<Integer> split;
				if (compList.size() < list.size())
					split = compList;
				else if (compList.size() > list.size())
					split = list;
				else {
					if (compList.get(0) < list.get(0))
						split = compList;
					else
						split = list;
				}

				// In trivial case
				if (split.size() == 1)
					continue;

				// Assume unique set of descendants on internal
				splitSet.add(split);
			}
		}
		return splitSet;
	}

	// For the sorting of HashMap with respect to its values
	public static LinkedHashMap<Integer, Integer> sortHashMapByValues(LinkedHashMap<Integer, Integer> passedMap) {

		LinkedHashMap<Integer, Integer> sortedMap = new LinkedHashMap<Integer, Integer>();
		List<Integer> mapValues = new ArrayList<Integer>(passedMap.values());
		Collections.sort(mapValues);
		Collections.reverse(mapValues);

		for (int i = 0; i < mapValues.size(); i++) {
			Integer val = mapValues.get(i);
			for (Integer key : passedMap.keySet()) {
				Integer mapValue = passedMap.get(key);
				if (val == mapValue) {
					passedMap.remove(key);
					sortedMap.put(key, val);
					break;
				}
			}
		}
		return sortedMap;
	}

	/*
	 * // This method will make frequency distribution over tree parameter public
	 * static LinkedHashMap<Integer,Integer> getTreeDist(ArrayList<Sample> samples,
	 * double burnin) throws IOException {
	 * 
	 * int startIndex = (int) (samples.size() * burnin);
	 * 
	 * // Leaf map, descendants Maps and frequency distribution over tree parameter
	 * HashMap<String, Integer> name2id = new HashMap<String, Integer>();
	 * HashMap<Integer, String> id2name = new HashMap<Integer, String>();
	 * LinkedHashMap<Integer,HashSet<List<Integer>> > desMaps = new
	 * LinkedHashMap<Integer,HashSet<List<Integer>>>();
	 * LinkedHashMap<Integer,Integer> treeDist = new
	 * LinkedHashMap<Integer,Integer>();
	 * 
	 * // Constructing trees and converting them into split (Descendants) maps
	 * for(int i=startIndex ; i<samples.size() ; i++) {
	 * 
	 * String treeStr = samples.get(i).tree; Tree tree =
	 * NewickReader.readNewickTreeStr( treeStr );
	 * 
	 * if(i==startIndex) setLeafMaps(tree.getLeafNames(true),name2id, id2name);
	 * 
	 * HashSet<List<Integer>> splitsSet = getDescendantsSet(tree,name2id,id2name);
	 * desMaps.put(i, splitsSet); }
	 * 
	 * 
	 * // Making frequency distribution over tree parameter int totalCases = 0;
	 * for(int i=startIndex ; i < samples.size() ; i++) {
	 * 
	 * totalCases++; HashSet<List<Integer>> desSet = desMaps.get(i);
	 * if(treeDist.size() == 0) { treeDist.put(i, 1); continue; }
	 * 
	 * // Matching with all trees in the distribution boolean found = false;
	 * for(Integer key:treeDist.keySet()) { HashSet<List<Integer>> binSet =
	 * desMaps.get(key); if( compareDescSet(desSet,binSet) ) { treeDist.put(key,
	 * treeDist.get(key) + 1); found = true; break; } } if(!found) treeDist.put(i,
	 * 1);
	 * 
	 * }
	 * 
	 * return treeDist;
	 * 
	 * }
	 */

	// This method will set the leave maps for the given tree
	public static void setLeafMaps(String[] leaves, HashMap<String, Integer> name2id,
			HashMap<Integer, String> id2name) {
		int id = 0;
		for (String leaf : leaves) {
			if (leaf != null) {
				name2id.put(leaf, id);
				id2name.put(id++, leaf);
			}
		}
	}

	// This method will compare the descendants maps induced by two trees for the
	// equality purposes
	public static boolean compareDescSet(HashSet<List<Integer>> desSet1, HashSet<List<Integer>> desSet2) {
		for (List<Integer> element : desSet1)
			if (!desSet2.contains(element))
				return false;
		return true;
	}

	// Parse tree object to label each node with the descendants below it.
	// The root and leaves are considered to be trivial cases
	// The descendants are returned in Hashable format
	public static HashSet<List<Integer>> getDescendantsSet(Tree tree, HashMap<String, Integer> name2id,
			HashMap<Integer, String> id2name) {

		HashSet<List<Integer>> descendantsSet = new HashSet<List<Integer>>();

		for (int i = 0; i < tree.nodes.size(); i++) {

			Node node = tree.nodes.get(i);
			if (node.nchildren == 0 || node.parent == null) {
				continue;
			} else {

				// Get descendants under this node
				ArrayList<String> descendants = new ArrayList<String>();
				getDescendants(node, descendants);

				// Get into the list of integers and sort them
				ArrayList<Integer> list = new ArrayList<Integer>();
				for (String id : descendants)
					list.add(name2id.get(id));
				Collections.sort(list);

				// Assume unique set of descendants on internal
				descendantsSet.add(list);
			}
		}
		return descendantsSet;
	}

	// Recursive method to get descendants below 'node'
	static void getDescendants(Node node, ArrayList<String> list) {
		for (int i = 0; i < node.nchildren; i++) {
			Node child = node.children.get(i);
			if (child.nchildren == 0) {
				list.add(child.name);
			} else
				getDescendants(child, list);
		}
	}

}