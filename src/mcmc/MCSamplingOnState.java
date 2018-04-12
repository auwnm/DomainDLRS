package mcmc;

import java.io.IOException;
import java.util.ArrayList;

import phylogeny.Reconciliation;
import rateModel.RateModel;
import seqEvolution.RKSeqEvoModel;
import seqEvolution.SubstitutionModel;
import Main.Parameters;
import Main.TuningParameters;
import common.Log;
import common.LogDouble;
import common.PRNG;
import dlModel.BirthDeathProbs;
import dlModel.DLModel;

public class MCSamplingOnState {

	static int NO_OF_SAMPLES = 100;

	static Reconciliation geneRecon;
	static Reconciliation[] domainRecon;

	static DLModel geneDLModel;
	static DLModel[] domDLModel;
	static RateModel[] rateModel;

	public static StateLikelihood compute(Parameters pi, PRNG prng, boolean useRootEdge, Log log) throws IOException {

		// Initializing variables for this state
		LogDouble product;
		ArrayList<LogDouble> productList = new ArrayList<LogDouble>();

		// Computing MPR reconciliations between domain,gene and species tree
		geneRecon = new Reconciliation(pi.geneTree, pi.speciesTree, pi.gene2species);
		domainRecon = new Reconciliation[pi.numberOfDomains];
		for (int i = 0; i < pi.numberOfDomains; i++)
			domainRecon[i] = new Reconciliation(pi.domainTree[i], pi.geneTree, pi.domain2genes[i]);

		// DL Model for gene and domain trees
		double[] speciesExtinctionProb = BirthDeathProbs.calcExtinctionProb(pi.speciesTree, pi.geneBDRates);
		geneDLModel = new DLModel(geneRecon, pi.geneBDRates, speciesExtinctionProb, useRootEdge, prng, log);
		domDLModel = new DLModel[pi.numberOfDomains];
		for (int r = 0; r < pi.numberOfDomains; r++)
			domDLModel[r] = new DLModel(domainRecon[r], pi.domBDRates[r], null, useRootEdge, prng, log);

		double[] extinctionProb;
		LogDouble domDLLikelihood, domRTLikelihood;
		LogDouble geneDLModelLikelihood = geneDLModel.birthDeathTreePrior();

		// Rate Model initialized
		RateModel[] rateModel = new RateModel[pi.numberOfDomains];
		for (int i = 0; i < pi.numberOfDomains; i++)
			rateModel[i] = new RateModel(pi.domEdgeRates[i], useRootEdge);

		for (int n = 0; n < NO_OF_SAMPLES; n++) {

			// Initializing Product variable
			product = new LogDouble(1.0);

			// Sample times on gene tree
			geneDLModel.sampleRealisation();

			// Calling DLR Model for each domain tree
			for (int r = 0; r < pi.numberOfDomains; r++) {

				extinctionProb = BirthDeathProbs.calcExtinctionProb(pi.geneTree, pi.domBDRates[r]);
				domDLModel[r].updateDLModel(extinctionProb);
				domDLLikelihood = domDLModel[r].birthDeathTreePrior();
				product.mult(domDLLikelihood);

				domDLModel[r].sampleRealisation();
				rateModel[r].update(domainRecon[r].guest);
				domRTLikelihood = rateModel[r].getLogLikelihood();
				product.mult(domRTLikelihood);

				productList.add(product);
			}

		}

		// Initializing statelikelihood vector
		StateLikelihood stateLikelihood = new StateLikelihood(pi.numberOfDomains);
		stateLikelihood.geneDLModel = geneDLModelLikelihood;
		stateLikelihood.domainDLRModel = LogDouble.logSumOfExp(productList).div(NO_OF_SAMPLES);

		if (TuningParameters.useNativeCode) {
			for (int r = 0; r < pi.numberOfDomains; r++)
				stateLikelihood.domainSeqEvoModel[r] = RKSeqEvoModel.getSeqLikelihood(pi.domainTree[r], pi.domMSA[r]);
		} else {
			for (int r = 0; r < pi.numberOfDomains; r++) {
				SubstitutionModel seqEvomodel = new SubstitutionModel(pi.domMSA[r], pi.Q, pi.domainTree[r],
						useRootEdge);
				stateLikelihood.domainSeqEvoModel[r] = seqEvomodel.getLogLikelihood();
			}
		}
		stateLikelihood.updateOverall();

		return stateLikelihood;
	}

	////////

	/*
	 * 
	 * 
	 * // Recording maximum value //LogDouble maxLikelihood = new LogDouble(0.0);
	 * //if(Product.greaterThanOrEquals(maxLikelihood)) // maxLikelihood = Product;
	 * // Recording products // productList.add(Product);
	 * 
	 * // Sort the products and summing them up. //LogDouble Sum =
	 * LogDouble.logSumOfExp(productList); // Divided by total number of iteration
	 * //Sum.div((double) NO_OF_SAMPLES); //stateLikelihood.mult(Sum);
	 * 
	 * 
	 * // To estimate the performance long startTime = System.nanoTime();
	 * 
	 * 
	 * 
	 * long stopTime = System.nanoTime(); System.out.println(
	 * "\nApprox. computational time (sec) for MC Sampling on state =  " + (double)
	 * (stopTime - startTime) / Math.pow(10,9));
	 * 
	 * // Initialize state likelihood LogDouble stateLikelihood = new LogDouble
	 * (1.0);
	 * 
	 * 
	 * 
	 */

	// Writing parameters
	// if(WriteStateParameters)
	// log.writeParameters(pi,geneRecon,domainRecon,true);
	// static boolean WriteStateParameters = true;

	// Log stateLog = new Log(mainSystem.CWD + "log.state.out");
	// Logging the probabilities for this iteration
	// stateLog.writeIteration(n,pi,geneDLModelLikelihood,domainDLModelLikelihood,rateModelLikelihoods,Product);
	// stateLog.close();

	// System.out.println("Maximum Likelihood = " + maxLikelihood);
	// System.out.println("Expected Likelihood = " + Sum.getLogValue());

	/*
	 * // Performing sequence evolution at domain Level LogDouble dataLikelihood =
	 * new LogDouble(1.0); SubstitutionModel smodel[] = new
	 * SubstitutionModel[pi.numberOfDomains]; for(int i=0; i < pi.numberOfDomains ;
	 * i++ ) { smodel[i] = new SubstitutionModel(pi.domMSA[i], pi.Q,
	 * pi.domainTree[i],false); dataLikelihood.mult(smodel[i].getLogLikelihood());
	 * System.out.println("Data likelihood for domain tree " + (i+1) + ":" +
	 * smodel[i].getLogLikelihood()); }
	 * 
	 * System.out.println("Data likelihood = " + dataLikelihood.getLogValue());
	 * 
	 * stateLikelihood.mult(dataLikelihood);
	 * 
	 */

	/*
	 * //rateModelLikelihoods[i] =
	 * RateModel.getLogLikelihood1(pi.domEdgeRates[i],pi.domainTree[i],log); static
	 * boolean FirstMCMCIteration = false; // If this is the first MCMC iteration,
	 * initialize the birthdeath & rate parameters. if(FirstMCMCIteration) {
	 * ParametersHandler.estimateBDRates(pi, geneRecon, domainRecon);
	 * ParametersHandler.estimateEdgeRates(pi,geneRecon,domainRecon,log);
	 * //ParametersHandler.InitializeParametersWithActualEstimates(pi);
	 * System.out.println("Initialisation done!"); }
	 */

}
