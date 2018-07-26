package mcmc;

import java.io.IOException;
import java.util.ArrayList;

import common.Log;
import common.LogDouble;
import common.PRNG;

import proposers.Proposer;
import proposers.ProposerSelector;
import proposers.ProposerStatistics;

import Main.OutputHandler;
import Main.Parameters;
import Main.TuningParameters;


public class MCMCManager {

	private Log log;
	private PRNG prng;
	private Parameters pi;
	private ProposerSelector ps;
	private RealParameterUniformPrior edgeRatePriors;
	private boolean useRootEdge;
			
	private int iteration;
	private double percent;
	public ProposerStatistics runStats;
	private OutputHandler output;
	public ProposalAcceptor pa;
	public boolean hillClimbingExhausted;
	public StringBuilder bestSample;

	
	public MCMCManager(Parameters pi,RealParameterUniformPrior edgeRateDistParameterPrior, ProposerSelector ps, ProposalAcceptor proposalAcceptor, OutputHandler output,PRNG prng,boolean useRootEdge, Log log) {
		this.ps = ps;
		this.pi = pi;
		this.edgeRatePriors = edgeRateDistParameterPrior;
		this.log = log;
		this.prng = prng;
		this.output = output;
		this.pa = proposalAcceptor;
		this.bestSample = new StringBuilder();
		this.useRootEdge = useRootEdge;
		
		this.iteration = 0;
		this.percent = TuningParameters.mcmcMaxIteration / 100.0 ;
		this.runStats  = new ProposerStatistics("Run Statistics");
		
	}

	public void run() throws IOException {
		
		ArrayList<Proposer> proposersList = new ArrayList<Proposer>();
		StateLikelihood oldStateLikelihood,proposedStateLikelihood,bestStateLikelihood;
		
		
		// Running MCMC chain or Hillclimbing search 
		if(TuningParameters.runType.equalsIgnoreCase("MCMC"))
      		System.out.println("Running MCMC chain ...");
      	else
      		System.out.println("Running Hillclimbing technique ...");
      	
				
		// Initial likelihood estimate
		if(TuningParameters.useMCSamplingOnState) {
			oldStateLikelihood  = MCSamplingOnState.compute(pi, prng, useRootEdge, log);
			bestStateLikelihood = oldStateLikelihood;
		} else {
			HillClimbingOnState.SetHillClimbingParameters(prng,edgeRatePriors, useRootEdge, log);
			oldStateLikelihood = HillClimbingOnState.compute(pi,iteration);
			bestStateLikelihood = oldStateLikelihood;
		}
			
      	
      	
		// Main loop
		do {

			// Clear the propoers list
			proposersList.clear();
		   			
			// Selection of proposers & perturbation of state
			ps.SelectProposers(proposersList);
			for (Proposer proposer : proposersList) 
				proposer.cacheAndPerturb();
			
						
			// Likelihood of the proposed state
			if(TuningParameters.useMCSamplingOnState)
				proposedStateLikelihood = MCSamplingOnState.compute(pi, prng, useRootEdge, log);
			else
				proposedStateLikelihood = HillClimbingOnState.compute(pi,iteration);
						
			
			// Do Accept or reject
			boolean doAccept = false;
			LogDouble newPosteriorDensity = proposedStateLikelihood.getUnnormalizedDensity();
			LogDouble oldPosteriorDensity = oldStateLikelihood.getUnnormalizedDensity();
			doAccept = this.pa.acceptProposedState(newPosteriorDensity,oldPosteriorDensity,proposersList);
						 
			if (doAccept) {
				
				runStats.increment(true, "" + proposersList.size() + " used proposers");
				for (Proposer proposer : proposersList) 
					proposer.clearCache();
				oldStateLikelihood = proposedStateLikelihood;
				
				if (bestStateLikelihood.unnormalizedDensity.lessThan(proposedStateLikelihood.unnormalizedDensity)) {
					bestStateLikelihood = proposedStateLikelihood;
					this.bestSample = new StringBuilder();
					this.bestSample.append(	output.getGeneParametersStr(iteration, pi, bestStateLikelihood) );
					this.bestSample.append("\n");
					for(int r=0 ; r<pi.numberOfDomains ; r++) {
						this.bestSample.append( output.getDomainParametersStr(iteration,pi, r, bestStateLikelihood, true, true) );
						this.bestSample.append( "\n");
					}	
				}
				
				
			} else {
				runStats.increment(false, "" + proposersList.size() + " used proposers");
				for (Proposer proposer : proposersList) {
					proposer.restoreCache();	
				}
			}
									
			// Sample state for chain as per thinning factor
			if( (iteration % TuningParameters.thinningFactor) == 0 ) {
				output.writeGeneParameters(iteration, pi, oldStateLikelihood);
				for(int r=0 ; r<pi.numberOfDomains ; r++) 
					output.writeDomainParameters(iteration,pi, r, oldStateLikelihood, true, true);
					
				if( iteration %  (percent * 10.0) == 0 && iteration != 0 )
					System.out.print( "\n" + String.format("%d",(int) (iteration/percent)) + "% completed out of maximum allowable iterations ... ");
			}
			
			if ( this.pa.hasExhausted() ) {
				log.write("\nOptimum reached: " + TuningParameters.hillClimbingMaxTries + " iterations without improvement.");
				log.newLine();
				break;
			}
						
			//Updating log file
			log.commit();

			
		} while (++iteration < TuningParameters.mcmcMaxIteration);
		
	}
		
}











//////////////////////////////////////////////////
//Rate Model Perturbation
//////////////////////////////////////////////////
/*
 			log.write("\n\nBefore Perturbation: ");
			log.write(pi.domEdgeRates[0][0] + "\t" + pi.domEdgeRates[0][1] + "\t" + pi.domEdgeRates[1][0] + "\t" + pi.domEdgeRates[1][1] );

			for(int r=0; r<pi.numberOfDomains ; r++) {
				for(int j=0 ; j<2 ; j++) {
					if(pi.domEdgeRates[r][j]>5.0)
						System.out.println("Error");

				}
			}

			log.write("\nAcceptance:");
			log.write( pi.domEdgeRates[0][0] + "\t" + pi.domEdgeRates[0][1] + "\t" + pi.domEdgeRates[1][0] + "\t" + pi.domEdgeRates[1][1] );

			log.write("\nRejected:");
			log.write( pi.domEdgeRates[0][0] + "\t" + pi.domEdgeRates[0][1] + "\t" + pi.domEdgeRates[1][0] + "\t" + pi.domEdgeRates[1][1] );
			log.commit();




			for(int r=0; r<pi.numberOfDomains ; r++) {
				for(int j=0 ; j<2 ; j++) {
					if(pi.domEdgeRates[r][j]>5.0)
						System.out.println("Error");

				}
			}
			
   
 */






/*
////////////////////////////////////////////////////
// Testing Gene Tree before and after perturbation
////////////////////////////////////////////////////
// import phylogeny.NewickWriter;
if(ps.level == Parameters.GENE_LEVEL && ! proposersList.get(0).getProposerName().equalsIgnoreCase("geneBDRates")) {
	System.out.println("\n=======================================");
	System.out.println( "Before perturbation at " + iteration + "\t" + NewickWriter.getNewickString(pi.geneTree, pi.geneTree.vt,true) );
						
}	

// In case of acceptance 
if(ps.level == Parameters.GENE_LEVEL && ! proposersList.get(0).getProposerName().equalsIgnoreCase("geneBDRates")) {
	System.out.println( "Acceptance at" + "\t" + iteration + "\t" + NewickWriter.getNewickString(pi.geneTree, pi.geneTree.vt,true) );
	System.out.println("=======================================");
	//System.out.println( "Acceptance at" + "\t" + iteration + "\t" + oldStateLikelihood.show() + "\t"  + "proposed_lk = " + proposedStateLikelihood.show());
}



// In case of rejection
if(ps.level == Parameters.GENE_LEVEL && ! proposersList.get(0).getProposerName().equalsIgnoreCase("geneBDRates")) {
	System.out.println( "Rejection" + "\t" + oldStateLikelihood.show() + "\t"  + "proposed_lk = " + proposedStateLikelihood.show());
}	

// After restoring				 
if(ps.level == Parameters.GENE_LEVEL && ! proposersList.get(0).getProposerName().equalsIgnoreCase("geneBDRates")) {
	System.out.println( "After restoring at     " + iteration + "\t" + NewickWriter.getNewickString(pi.geneTree, pi.geneTree.vt,true) );
	System.out.println("\n=======================================");
}

*/











//output.writeGeneTopology(iteration, pi);
//output.writeDomainTopology(iteration, pi, r);

/*
if(ps.level == Parameters.GENE_LEVEL && ! proposersList.get(0).getProposerName().equalsIgnoreCase("geneBDRates")) {
	//System.out.println(NewickWriter.getNewickString(pi.geneTree, pi.geneTree.vt,true));
	System.out.print(oldStateLikelihood.show());
	System.out.println( "\t proposed_lk = " + proposedStateLikelihood.show()  + "\t" + "Rejection");
	//System.out.println(NewickWriter.getNewickString(pi.geneTree, pi.geneTree.vt,true) + "\t" + " *Rejection" + "\t" + iteration);
	//System.out.print(" *Rejection ");
}
*/

/*
if(ps.level == Parameters.GENE_LEVEL && ! proposersList.get(0).getProposerName().equalsIgnoreCase("geneBDRates")) {
	//System.out.println(NewickWriter.getNewickString(pi.geneTree, pi.geneTree.vt,true));
	//System.out.print(oldStateLikelihood.show());
	//System.out.println( "\t proposed_lk = " + proposedStateLikelihood.show() + "\t" + "Acceptence");
	System.out.println(NewickWriter.getNewickString(pi.geneTree, pi.geneTree.vt,true) + "\t" + " *Acceptence" + "\t" + iteration);
	
}
*/


//System.out.println(NewickWriter.getNewickString(pi.geneTree, pi.geneTree.vt,true) + "\t" + "Time of root node =" +  pi.geneTree.vt[pi.geneTree.root.id]);
//if(ps.level == Parameters.GENE_LEVEL && ! proposersList.get(0).getProposerName().equalsIgnoreCase("geneBDRates")) {
//	System.out.println(NewickWriter.getNewickString(pi.geneTree, pi.geneTree.vt,true));
//}	



/*
LogDouble a = proposedStateLikelihood.overall.divToNew(oldStateLikelihood.overall); 
for (Proposer proposer : proposersList) {
	if (!proposer.hasValidProposal()) {
		doAccept = false;
	}
	a.mult( proposer.getDensityRatio() );
}
// Decision of Acceptance or Rejection
LogDouble randomVal = new LogDouble(prng.nextDouble()); 
if( a.greaterThanOrEquals(randomVal)  )
	doAccept = true;
 */

/*
//double [] length1 = Arrays.copyOf(pi.domainTree[0].bl,pi.domainTree[0].bl.length); 
//double [] length2 = Arrays.copyOf(pi.domainTree[1].bl,pi.domainTree[1].bl.length); 	

if(!doAccept) {

	//pi.domainTree[0].bl = length1;
	//pi.domainTree[1].bl = length2;

	for(int i=0 ; i<pi.domainTree[0].bl.length ; i++)
		if(pi.domainTree[0].bl[i] != length1[i])
			System.out.println("Error1");

	for(int i=0 ; i<pi.domainTree[1].bl.length ; i++)
		if(pi.domainTree[1].bl[i] != length2[i])
			System.out.println("Error2");

}
 */

//System.out.print(oldStateLikelihood.show());
//System.out.println( "\t proposed_lk = " + proposedStateLikelihood.show() );


//pi.showState(0);
//pi.showState(0);
//oldStateLikelihood = MCSamplingOnState.compute(pi,log);
//proposedStateLikelihood = MCSamplingOnState.compute(pi,log);
//System.out.println(iteration + "\t ratio = " + a.getLogValue() +"\t rv = " + randomVal.getLogValue());


//double totime = 0.0;
//long startTime, stopTime ;
//log.newLine();
//log.write("Approx. mean computational time for single state in (sec) = " + ( totime / (double) TuningParameters.mcmcMaxIteration ) );
//startTime = System.nanoTime();
//stopTime = System.nanoTime();
//totime += (double) (stopTime - startTime) / Math.pow(10,9) ; 
