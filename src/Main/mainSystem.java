package Main;

import java.util.ArrayList;

import mcmc.HillClimbingAcceptor;
import mcmc.MCMCManager;
import mcmc.MetropolisHastingsAcceptor;
import mcmc.ProposalAcceptor;
import mcmc.RealParameterUniformPrior;

import proposers.LengthProposer;
import proposers.Proposer;
import proposers.ProposerSelector;
import proposers.RateProposer;
import proposers.TopologyProposer;

import common.Log;
import common.PRNG;
import common.RealInterval;

public class mainSystem {

	public static String CWD;
	public static int numberOfDomains;
	public static long startTime, stopTime;

	public static void main(String[] args) {

		try {

			// Parsing input parameters & note start time
			startTime = System.nanoTime();
			Input input = InputHandler.readInputArguments(args);
			if (input == null)
				return;
			CWD = input.path;
			numberOfDomains = input.numberOfDomains;

			// Initialize System logging and writing pre-run info
			Log log = new Log(CWD + input.outPrefix + ".inf");
			log.writePreRunInfo(input);

			// Loading native code libaray
			if (TuningParameters.useNativeCode)
				System.load("/Users/auwnm/Documents/cworkspace/RasmussenSeqEvo/Debug/libRasmussenSeqEvo.dylib");

			// Initialising pseudo-random number generator
			PRNG prng = ParametersHandler.getPRNG(TuningParameters.seed);

			// Setting Priors, We only have them for parameters which might cause issues.
			RealInterval priorRange = new RealInterval(0, 5.0, false, false);
			RealParameterUniformPrior edgeRatePriors = new RealParameterUniformPrior(priorRange, false);

			// Initialising parameters
			Parameters pi = new Parameters(numberOfDomains);
			ParametersHandler.initSubstitutionModel(pi, input.substitutionModel);
			ParametersHandler.readGeneMappingFile(pi, CWD, input.geneMapFile);
			ParametersHandler.readDomainMappingFiles(pi, CWD, input.domMappingFiles);
			ParametersHandler.readAlignmentFiles(pi, CWD, input.domAlignmentFiles);
			ParametersHandler.readSpeciesTree(pi, CWD, input.speciesTreeFile);
			ParametersHandler.initDomainTrees(pi, prng, TuningParameters.initNJDomainTree);
			ParametersHandler.initGeneTree(pi, input.geneTree, prng, TuningParameters.initNJGeneTree);
			ParametersHandler.initRates(pi, edgeRatePriors, prng, TuningParameters.useRootEdge, log);

			// Initialising Output Handler
			OutputHandler output = new OutputHandler(numberOfDomains, CWD, input.outPrefix);

			// Set domain perturbation weights
			if (TuningParameters.domainSelectorWeights == null)
				TuningParameters.setDomainPerturbationWeights(pi, false);

			// Initialising Proposer Selector
			ProposerSelector ps = new ProposerSelector(numberOfDomains, prng);

			// Gene Level Proposers
			ps.geneLevelProposers.add(new RateProposer(Parameters.BD_RATES, pi, Parameters.GENE_LEVEL, -1,
					TuningParameters.BDRatesCV, prng));
			ps.geneLevelProposers.add(
					new TopologyProposer(pi, Parameters.GENE_LEVEL, -1, TuningParameters.geneTreeMoveWeights, prng));

			// Domain Level Proposers
			for (int r = 0; r < numberOfDomains; r++) {
				ArrayList<Proposer> proposers = new ArrayList<Proposer>();
				proposers.add(new RateProposer(Parameters.BD_RATES, pi, Parameters.DOMAIN_LEVEL, r,
						TuningParameters.BDRatesCV, prng));
				proposers.add(new RateProposer(Parameters.EDGE_RATES, pi, Parameters.DOMAIN_LEVEL, r,
						TuningParameters.EdgeRatesCV, prng));
				proposers.add(
						new LengthProposer(pi, r, TuningParameters.branchLengthCV, prng, TuningParameters.useRootEdge));
				proposers.add(new TopologyProposer(pi, Parameters.DOMAIN_LEVEL, r, TuningParameters.domTreeMoveWeights,
						prng));
				ps.domLevelProposers.put(r, proposers);
			}

			// Enabling or Disabling the proposers
			ps.geneLevelProposers.get(0).setEnabled(true);

			if (TuningParameters.fixedGeneTree && input.geneTree != null)
				ps.geneLevelProposers.get(1).setEnabled(false);
			else
				ps.geneLevelProposers.get(1).setEnabled(true);

			for (int r = 0; r < numberOfDomains; r++) {
				ps.domLevelProposers.get(r).get(Parameters.BD_RATES).setEnabled(true);
				ps.domLevelProposers.get(r).get(Parameters.EDGE_RATES).setEnabled(true);
				ps.domLevelProposers.get(r).get(Parameters.LENGTHS).setEnabled(true);
				ps.domLevelProposers.get(r).get(Parameters.TOPOLOGY).setEnabled(true);
			}

			// Proposal Acceptor
			ProposalAcceptor pa;
			if (TuningParameters.runType.equalsIgnoreCase("HILLCLIMBING"))
				pa = new HillClimbingAcceptor(TuningParameters.hillClimbingMaxTries);
			else
				pa = new MetropolisHastingsAcceptor(prng);

			// Running MCMC chain on top of that
			MCMCManager mcmcManager = new MCMCManager(pi, edgeRatePriors, ps, pa, output, prng,
					TuningParameters.useRootEdge, log);
			mcmcManager.run();

			// Writing post run info & closing output & system log files
			log.writePostRunInfo(ps, mcmcManager.runStats);
			stopTime = System.nanoTime();
			log.write(String.format("Approx. computational time (sec) for %d = %f", TuningParameters.mcmcMaxIteration,
					(stopTime - startTime) / Math.pow(10, 9)));
			log.writeBestState(mcmcManager.bestSample);
			log.close();
			output.close();
			System.out.print("\ndone.");

		} catch (Exception e) {
			System.out.println(e.toString());
			e.printStackTrace();
		}
	}
}

/*
 * // Files & input variables String speciesTreeFile; String geneMapFile; int
 * numberOfDomains = 2; String substitutionModel = "JC69"; String []
 * domAlignmentFiles = new String[numberOfDomains]; String [] domMappingFiles =
 * new String[numberOfDomains];
 * 
 * 
 * // Specifying input speciesTreeFile = CWD + "hostTree.nwk"; geneMapFile = CWD
 * + "geneTree.map"; domMappingFiles[0] = CWD + "domainTree1.map";
 * domMappingFiles[1] = CWD + "domainTree2.map"; domAlignmentFiles[0] = CWD +
 * "seq_relaxed_domainTree1.fasta"; domAlignmentFiles[1] = CWD +
 * "seq_relaxed_domainTree2.fasta";
 */

// String lk = HillClimbingOnState.compute(pi, log).show();
// System.out.println(lk);

/*
 * 
 * HillClimbingOnState.compute(pi, log).show(); System.out.println();
 * ps.domLevelProposers.get(0).get(3).cacheAndPerturb();
 * System.out.println("name of proposer =" +
 * ps.domLevelProposers.get(0).get(3).getProposerName()); System.out.println();
 * HillClimbingOnState.compute(pi, log).show();
 * 
 */

// StateLikelihood sl1 = HillClimbingOnState.compute(pi, log);
// System.out.println("Hill Climbing answer = " + sl1.overallLikelihood());
// StateLikelihood sl2 = MCSamplingOnState.compute(pi, log);
// System.out.println("MC Sampling answer = " + sl2.overallLikelihood());
