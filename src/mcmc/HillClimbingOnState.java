package mcmc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import phylogeny.Node;
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
import dlModel.Realization;


public class HillClimbingOnState {

	static PRNG prng;
	static Log log;
	static boolean useRootEdge;
	
	static int MAX_HILL_TRIES = 5;
	static int MAX_NO_SAMPLES = 20;
	
	static double RESOLUTION   = 1e3; 
	static double ACCELERATION = 1.5;
	
	
	static RealParameterUniformPrior edgeRatePriors;
	
	static Reconciliation geneRecon;
	static Reconciliation [] domainRecon;

	static DLModel geneDLModel;
	static DLModel [] domDLModel;
	static RateModel [] rateModel;
		
		
	public static void SetHillClimbingParameters(PRNG prng,RealParameterUniformPrior edgeRatePriors,boolean useRootEdge,Log log) {
		HillClimbingOnState.edgeRatePriors = edgeRatePriors;
		HillClimbingOnState.prng = prng;
		HillClimbingOnState.log = log;
		HillClimbingOnState.useRootEdge = useRootEdge;		
	}
		
	
	public static StateLikelihood compute(Parameters pi,int iteration) throws IOException {
		
		// Computing MPR reconciliations between domain,gene and species tree
		geneRecon   = new Reconciliation(pi.geneTree,pi.speciesTree,pi.gene2species);
		domainRecon = new Reconciliation[pi.numberOfDomains];
		for(int r=0 ; r < pi.numberOfDomains ; r++) 
			domainRecon[r] = new Reconciliation(pi.domainTree[r],pi.geneTree, pi.domain2genes[r]);
		
		// Gene DL Model 
		double [] extinctionProb = BirthDeathProbs.calcExtinctionProb(pi.speciesTree,pi.geneBDRates);
		geneDLModel = new DLModel(geneRecon,pi.geneBDRates,extinctionProb,useRootEdge,prng,log);
		LogDouble geneDLModelLikelihood = geneDLModel.birthDeathTreePrior();
		
		// Setting up domain DL Model
		domDLModel = new DLModel [pi.numberOfDomains];
		for(int r=0 ; r < pi.numberOfDomains ; r++)
			domDLModel[r]  = new DLModel(domainRecon[r], pi.domBDRates[r],null,useRootEdge,prng,log);
		
		// Setting up domain Rate Model 
		rateModel = new RateModel[pi.numberOfDomains];
		LogDouble meanPrior,cvPrior;
		for(int r=0; r<pi.numberOfDomains ; r++) {
			meanPrior = edgeRatePriors.getPriorProbability(pi.domEdgeRates[r][0]);
			  cvPrior = edgeRatePriors.getPriorProbability(pi.domEdgeRates[r][1]);
			if( meanPrior.isZero() || cvPrior.isZero() ) {
				StateLikelihood stateLikelihood = new StateLikelihood(pi.numberOfDomains);
				stateLikelihood.unnormalizedDensity = new LogDouble(0.0);
				return stateLikelihood;
			}
			rateModel[r] = new RateModel(pi.domEdgeRates[r],useRootEdge); 
		}
	
		// Sampling step
		LogDouble newLikelihood;
		LogDouble currentLikelihood = new LogDouble(0.0);
		LogDouble domDLLikelihood,domRTLikelihood;
		
		// Realization mechanism
		Realization [] domainReal = new Realization[pi.numberOfDomains];
		Realization geneReal = new Realization();
		for(int r=0 ; r < pi.numberOfDomains ; r++) 
			domainReal[r] = new Realization();
				
		
		boolean correctlySampled = true;
		for (int n=0 ; n< MAX_NO_SAMPLES ; n++ ) {

			// Start sampling
			correctlySampled = true;
			newLikelihood = new LogDouble(1.0);
			correctlySampled &= geneDLModel.sampleRealisation();
								
			// Calling DLR Model for each domain tree
			for(int r=0 ; r < pi.numberOfDomains ; r++) {
								
				extinctionProb = BirthDeathProbs.calcExtinctionProb(pi.geneTree,pi.domBDRates[r]);
				domDLModel[r].updateDLModel(extinctionProb);
				domDLLikelihood = domDLModel[r].birthDeathTreePrior();
				newLikelihood.mult(domDLLikelihood);
							
				correctlySampled &= domDLModel[r].sampleRealisation();
				rateModel[r].update(domainRecon[r].guest);
				domRTLikelihood = rateModel[r].getLogLikelihood();
				newLikelihood.mult(domRTLikelihood);
			}
						
			if(!correctlySampled) 
				continue;
			
			if(newLikelihood.greaterThan(currentLikelihood)) {
				geneReal.Cache(geneRecon);
				for(int r=0 ; r < pi.numberOfDomains ; r++)
					domainReal[r].Cache(domainRecon[r]);
				currentLikelihood = newLikelihood;
			} 
			
		}
			
		if(!correctlySampled) {
			StateLikelihood stateLikelihood = new StateLikelihood(pi.numberOfDomains);
			stateLikelihood.unnormalizedDensity = new LogDouble(0.0);
			return stateLikelihood;
		} else {
			geneReal.RestoreCache(geneRecon);
			for(int r=0 ; r < pi.numberOfDomains ; r++) 
				domainReal[r].RestoreCache(domainRecon[r]);
		}
				
		// Start hill climbing on that
		// Moveable nodes for hill climbing technique
		LinkedHashMap<Integer,ArrayList<Node>> domDups  = getDomainDuplications(domainRecon);
		ArrayList<Node> genDups = getGeneDuplications(geneRecon);
		
		// Starting likelihood business
		int iterations=0;
		currentLikelihood = new LogDouble(0.0);

		// Outermost loop
		do {
			
			// Maximize likelihood  
			newLikelihood = maximizeLikelihood(pi,domDups,genDups,log);
			if(newLikelihood.greaterThan(currentLikelihood)) 
				currentLikelihood = newLikelihood;
			else 				
				break;
		
		} while( ++iterations < MAX_HILL_TRIES);
		
		
		// Initializing statelikelihood vector
		StateLikelihood stateLikelihood = new StateLikelihood(pi.numberOfDomains);
		stateLikelihood.geneDLModel = geneDLModelLikelihood;
		stateLikelihood.domainDLRModel = currentLikelihood;
		
		//Performing sequence evolution
		if(TuningParameters.useNativeCode) {		
			for(int r=0 ; r < pi.numberOfDomains ; r++) 
				stateLikelihood.domainSeqEvoModel[r] = RKSeqEvoModel.getSeqLikelihood(pi.domainTree[r], pi.domMSA[r]);
		} else {
			for(int r=0; r < pi.numberOfDomains ; r++ ) {
				SubstitutionModel seqEvomodel = new SubstitutionModel(pi.domMSA[r], pi.Q, pi.domainTree[r],useRootEdge);
				stateLikelihood.domainSeqEvoModel[r] = seqEvomodel.getLogLikelihood();
			}
		}
		stateLikelihood.updateOverall();
		
		
		return stateLikelihood;
	}

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/*
	// Debuging code:
	//pi.domainTree[0].vt = null;
	//pi.domainTree[1].vt = null;
	
	//domainRecon[0].host.show(false);
	//System.out.println();
	//domainRecon[0].guestPrime.show(false);
	//System.out.println();
	//domainRecon[0].show();
	
	// Debug rate model
		log.newLine();
		log.write(iteration +"\t" + "Gamma( " + rateModel[0].getMean() + " , " + rateModel[0].getCV() + " )"  + "\t");
		//log.newLine();
		for(int i=0 ; i<pi.domainTree[0].nnodes ; i++) {
			if(rateModel[0].logLikelihoods[i] != null) {
				//log.write(rateModel[0].length[i] + "\t\t" + rateModel[0].time[i]  + "\t\t" + (rateModel[0].length[i] / rateModel[0].time[i])  + "\t Loglikelihood = " + rateModel[0].logLikelihoods[i].getLogValue());
				log.write((rateModel[0].length[i] / rateModel[0].time[i]) + "\t" );
				//log.newLine();
			}
		}
	
		
	// Checking for valid realization
	geneDLModel.isValidRealization(geneRecon);
	for(int r=0 ; r < pi.numberOfDomains ; r++)
		domDLModel[r].isValidRealization(domainRecon[r]);
	
	pi.showState();
	pi.geneTree.show(false);
	for(int k=0 ; k<100 ; k++) {
	 geneDLModel.sampleRealisation();
	 System.out.println(pi.geneTree.vt[8]);
	}
	
	long startTime = System.nanoTime();
	seqEvomodel = new SubstitutionModel[pi.numberOfDomains];
	for(int r=0; r < pi.numberOfDomains ; r++ ) {
		seqEvomodel[r] = new SubstitutionModel(pi.domMSA[r], pi.Q, pi.domainTree[r],false);
		stateLikelihood.domainSeqEvoModel[r] = seqEvomodel[r].getLogLikelihood();
		dataLikelihood.mult(seqEvomodel[r].getLogLikelihood());
	}

	long stopTime = System.nanoTime();
	System.out.println( "\nApprox. computational time (sec) = " +  (double) (stopTime - startTime) / Math.pow(10,9));
	
	for(int r=0 ; r < pi.numberOfDomains ; r++) {
		assert(domDLModel[r].isValidRealization(domainRecon[r]));
	}
	*/
	
	
	// Calculating complete likelihood score when time on gene tree has changed.
	public static LogDouble getCompleteLikelihood(Parameters pi, Log log) throws IOException {

		double [] extinctionProb;
		LogDouble DLModelLikelihood, RateModelLikelihood;
		LogDouble completeLikelihood = new LogDouble (1.0);

		// DLR Model Likelihood for each domain tree 
		for(int r=0 ; r < pi.numberOfDomains ; r++) {
			extinctionProb     =  BirthDeathProbs.calcExtinctionProb(pi.geneTree,pi.domBDRates[r]);
			domDLModel[r].updateDLModel(extinctionProb);
			DLModelLikelihood	= domDLModel[r].birthDeathTreePrior();
			RateModelLikelihood	= rateModel[r].getLogLikelihood();

			completeLikelihood.mult(DLModelLikelihood);
			completeLikelihood.mult(RateModelLikelihood);
		}

		return completeLikelihood;
	}


	public static LogDouble maximizeLikelihood(Parameters pi ,LinkedHashMap<Integer,ArrayList<Node>> domDups,ArrayList<Node> genDups, Log log) throws IOException {

		ArrayList<Node> nodes;
		double stepSize = geneRecon.host.getPeakTime() / RESOLUTION ;


		for(int r=0 ; r<pi.numberOfDomains ; r++) {

			if(!domDups.containsKey(r))
				continue;

			nodes = domDups.get(r);
			for(Node node: nodes) 
				maximizeDomainNode(node,r,stepSize);

		}

		for(int i=0; i<genDups.size() ; i++) 
			maximizeGeneNode(genDups.get(i),pi,stepSize,log);	


		return getCompleteLikelihood(pi,log);
	}



	// This method will maximize the likelihood score related to domain node
	public static void maximizeDomainNode(Node node, int r, double delta) {

		boolean moved,recovered;
		LogDouble bestLikelihood,currentLikelihood,tempLikelihood;
		double [] candidate = {-ACCELERATION , ACCELERATION };

		do {
			// Likelihood of the current state 
			currentLikelihood  = rateModel[r].getLogLikelihood();

			// Picking the best possible state
			bestLikelihood = new LogDouble(currentLikelihood);
			for(int c=0 ; c<2 ; c++) {
				moved = moveDomainNode(node, delta * candidate[c], domainRecon[r],rateModel[r]);
				if(moved) {
					tempLikelihood = rateModel[r].getLogLikelihood();
					if(tempLikelihood.greaterThan(bestLikelihood)) {
						bestLikelihood  = tempLikelihood;
						delta *= candidate[c]; 
						if(c==1) {
							candidate[0] *= -1.0 ;
							candidate[1] *= -1.0 ;
						}
						break;
					}	
					recovered = moveDomainNode(node, -1.0 * delta * candidate[c], domainRecon[r],rateModel[r]);
					assert(recovered);
				}
			}

		}while(bestLikelihood.greaterThan(currentLikelihood));

	}

	// This method will maximize the likelihood score related to gene node and all associated domain nodes
	public static void maximizeGeneNode(Node gnode, Parameters pi, double delta, Log log) throws IOException {

		boolean moved,recovered;
		LogDouble bestLikelihood,currentLikelihood,tempLikelihood;
		double [] candidate = {-ACCELERATION , ACCELERATION };

		do {

			// Likelihood of the current state 
			currentLikelihood = getCompleteLikelihood(pi,log);

			// Picking the best possible state
			bestLikelihood = new LogDouble(currentLikelihood);
			for(int c=0 ; c<2 ; c++) {
				moved = moveGeneNode(gnode,delta * candidate[c],geneRecon, domainRecon,rateModel);
				if(moved) {
					tempLikelihood  = getCompleteLikelihood(pi,log);
					if(tempLikelihood.greaterThan(bestLikelihood)) {
						bestLikelihood = tempLikelihood;
						delta *= candidate[c];
						if(c==1) {
							candidate[0] *= -1.0 ;
							candidate[1] *= -1.0 ;
						}
						break;
					}	
					recovered =  moveGeneNode(gnode, -1.0 * delta * candidate[c],geneRecon, domainRecon,rateModel); 
					assert(recovered);
				}
			}

		}while(bestLikelihood.greaterThan(currentLikelihood));

	}

	// This method will move gene node with all its reconciled associated domain nodes 
	// TODO: Optimize this function only for moveable nodes ....
	private static boolean moveGeneNode(Node gnode,double delta,Reconciliation geneRecon,Reconciliation [] domainRecon,RateModel [] rateModel ) {

		double t;
		double [] bounds;
		double [] validBounds = getTimeBoundry(gnode,geneRecon);

		ArrayList<Node> nodes; 
		LinkedHashMap<Integer,ArrayList<Node>> influenced = new LinkedHashMap<Integer,ArrayList<Node>>();
		for(int r=0 ; r < domainRecon.length ; r++ ) {
			nodes = new ArrayList<Node>();
			for(Node dnode: domainRecon[r].guestPrime.nodes) {			
				if( domainRecon[r].host2guestPrime[gnode.id][dnode.id] == 1 && domainRecon[r].eventMapPrime[dnode.id] == Reconciliation.EVENT_SPEC) {

					bounds = getTimeBoundry(dnode,domainRecon[r]); 

					if(bounds[0] > validBounds[0])
						validBounds[0] = bounds[0];
					if(bounds[1] < validBounds[1])
						validBounds[1] = bounds[1];

					nodes.add(dnode);
				}
			}
			influenced.put(r, nodes);	
		}

		// Checking if this move is valid or not
		boolean isValidMove = true;
		for(int r=0 ; r < domainRecon.length ; r++ ) {
			if(!influenced.containsKey(r)) 
				continue;
			for (Node dnode: influenced.get(r)) {
				t = domainRecon[r].guestPrime.vt[dnode.id];
				if( t+delta >= validBounds[1]  || t+delta <= validBounds[0] ) 
					isValidMove = false;
			}
		}
		t  = geneRecon.guestPrime.vt[gnode.id];
		if( t+delta >= validBounds[1]  || t+delta <= validBounds[0] ) 
			isValidMove = false;


		// If it is valid move then commit it
		if(isValidMove) {
			for(int r=0 ; r < domainRecon.length ; r++ ) {
				if(!influenced.containsKey(r)) 
					continue;
				for (Node dnode: influenced.get(r)) {
					t = domainRecon[r].guestPrime.vt[dnode.id];
					domainRecon[r].guestPrime.vt[dnode.id] = t+delta ;
					if(!dnode.isImplied) {
						Node tempNode = domainRecon[r].guest.nodes.get(dnode.id);
						domainRecon[r].guest.vt[tempNode.id] = t+delta;
						rateModel[r].update(domainRecon[r].guest);
					}

				}
			}
			geneRecon.guestPrime.vt[gnode.id] = t+delta ;
			if(!gnode.isImplied)
				geneRecon.guest.vt[gnode.id] = t+delta;
		} else
			return false;

		return true;
	}

	// This method will move only single domain node
	private static boolean moveDomainNode(Node node,double delta,Reconciliation domainRecon,RateModel rateModel ) {

		double [] validBounds = getTimeBoundry(node,domainRecon);

		// Checking if it move is doable
		double t  = domainRecon.guestPrime.vt[node.id];
		if( t+delta >= validBounds[1]  || t+delta <= validBounds[0] ) 
			return false;


		// Updating in prime tree i.e. tree with implied nodes
		// assuming same duplication id in both the trees
		// domainRecon.guestPrime.vt[node.id] = t+delta;
		domainRecon.guestPrime.vt[node.id] = t+delta ;


		// update rate mode
		if(!node.isImplied) {
			Node tempNode = domainRecon.guest.nodes.get(node.id);
			domainRecon.guest.vt[tempNode.id] = t+delta;
			rateModel.update(domainRecon.guest);
		}

		return true;
	}

	// This method will give the time boundary of a node
	private static double [] getTimeBoundry(Node node,Reconciliation recon) {

		double [] bounds = new double [2];

		double upperTime  = node.isRoot()?recon.guestPrime.getPeakTime():recon.guestPrime.vt[node.parent.id];

		double lowerTime = 0.0;
		for(int i=0 ; i<node.nchildren ; i++) {
			double childTime = recon.guestPrime.vt[node.children.get(i).id];
			if( childTime > lowerTime )
				lowerTime = childTime;
		}

		bounds[0] = lowerTime;
		bounds[1] = upperTime;
		return bounds;
	}


	// This method will extract the domain duplication nodes.
	public static LinkedHashMap<Integer,ArrayList<Node>> getDomainDuplications(Reconciliation [] domainRecon) {
		ArrayList<Node> nodes;
		LinkedHashMap<Integer,ArrayList<Node>> dups= new LinkedHashMap<Integer,ArrayList<Node>>();
		for(int r=0; r<domainRecon.length ; r++) {
			nodes = new ArrayList<Node>();
			for(Node dnode: domainRecon[r].guestPrime.nodes) {
				if( domainRecon[r].eventMapPrime[dnode.id]  == Reconciliation.EVENT_DUPL  ) 
					nodes.add(dnode);
			}
			dups.put(r, nodes);
		}
		return dups;
	}


	// This method will extract the gene duplication nodes.
	public static ArrayList<Node> getGeneDuplications(Reconciliation geneRecon) {

		ArrayList<Node> dups= new ArrayList<Node>();
		for(Node gnode: geneRecon.guestPrime.nodes) {
			if( geneRecon.eventMapPrime[gnode.id]  == Reconciliation.EVENT_DUPL  ) 
				dups.add(gnode);
		}
		return dups;
	}
}
