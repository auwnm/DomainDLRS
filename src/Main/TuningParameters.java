package Main;

public class TuningParameters {

	public static String seed = null; // PRNG seed. Default: Random seed
	public static String runType = "MCMC"; // Type of run. Valid values are MCMC and HILLCLIMBING. (default = "MCMC")
	public static int thinningFactor = 100; // Thinning factor, i.e., sample output every n-th iteration (default = 100)
	public static long mcmcMaxIteration = 1000000; // Maximum number of iterations for MCMC chain (default = 1 million)
	public static long hillClimbingMaxTries = 100; // Maximum number of tries for HillClimbing (default = 100)
	public static double minBranchLength = 1E-3; // Minimum branch length, branch length approaching to zero causes rate
													// model unstability

	public static boolean useNativeCode = false; // Usage of native code for speedup purpose (default = false)
	public static boolean useRootEdge = true; // If true, utilises the root arc ("stem") branch length when computing
												// model likelihood; if false, discards the root arc. (default = false)
	public static boolean useMCSamplingOnState = false; // Usage of Monte Carlo Sampling on mcmc state for likelihood
														// computation (default = false)

	public static boolean initNJGeneTree = true; // Constructing initial gene tree by using NJ algorithm on gene random
													// msa (default = true)
	public static boolean initNJDomainTree = true; // Constructing initial domain trees with their lengths by using NJ
													// algorithm (default = true)
	public static boolean fixedGeneTree = false; // Fixing the gene tree parameter (default = false)

	public static double BDRatesCV = 0.25; // Proposal distribution variance of birth death rate parameters.
	public static double EdgeRatesCV = 0.25; // Proposal distribution variance of Edge rate parameters.
	public static double branchLengthCV = 0.30; // Proposal distribution variance of branch lengths parameters.
	public static double[] geneTreeMoveWeights = { 0.45, 0.30, 0.20, 0.05 }; // Branch swap operation is carried out as
																				// [NNI,RNNI,SPR,Rerooting].
	public static double[] domTreeMoveWeights = { 0.35, 0.30, 0.30, 0.05 }; // Branch swap operation is carried out as
																			// [NNI,RNNI,SPR,Rerooting].

	public static double domainLevelWeight = 0.80; // Governs how often domain level parameters perturbation will be
													// performed.
	public static double[] domainSelectorWeights = { 0.30, 0.70 }; // Governs how often a particular domain is selected
																	// for perturbation. //{0.20,0.20,0.60};

	public static double[] geneParametersWeights = { 0.55, 0.45 }; // Governs how often a particular gene parameter
																	// [BDRates,Topology] will be selected for
																	// perturbation.
	public static double[] domParametersWeights = { 0.15, 0.15, 0.50, 0.20 }; // Governs how often a particular domain
																				// parameter
																				// [BDRates,EdgeRates,Lengths,Topology]
																				// will be selected for perturbation.
	public static double[] noProposersWeights = { 0.30, 0.50, 0.15, 0.05 }; // Governs how often 1,2,3 and all 4
																			// simultaneous proposers will be activated
																			// for performing a state change. (No more
																			// than 4 may be specified)
	public static double[] lengthsSelectorWeights = { 0.40, 0.35, 0.15, 0.10 }; // Governs how often 1,2,... branch
																				// lengths will be perturbed
																				// simultaneously, e.g., [0.5,0.5] for
																				// an equal chance of 1 or 2 branch
																				// lengths.

	// Setter method
	public static void setDomainPerturbationWeights(Parameters pi, boolean uniform) {

		int numberOfDomains = pi.numberOfDomains;
		domainSelectorWeights = new double[numberOfDomains];

		if (numberOfDomains == 1) {
			domainSelectorWeights[0] = 1.0;
			return;
		}

		if (uniform) {
			for (int i = 0; i < numberOfDomains; i++)
				domainSelectorWeights[i] = 1.0 / numberOfDomains;
			return;
		}

		int totalNodes = 0;
		for (int i = 0; i < numberOfDomains; i++)
			totalNodes += pi.domainTree[i].nnodes;

		for (int i = 0; i < numberOfDomains; i++)
			domainSelectorWeights[i] = (double) pi.domainTree[i].nnodes / (double) totalNodes;
	}

}
