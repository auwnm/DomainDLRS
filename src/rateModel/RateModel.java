package rateModel;

import common.LogDouble;
import phylogeny.Node;
import phylogeny.Tree;

public class RateModel {

	private double mean;
	private double cv;
	private GammaDistribution edgeRatePD;
	private boolean useRootEdge;

	public double[] length;
	public double[] time;
	public LogDouble[] logLikelihoods;

	// Constructor
	public RateModel(double[] edgeRates, boolean useRootEdge) {
		this.mean = edgeRates[0];
		this.cv = edgeRates[1];
		this.edgeRatePD = new GammaDistribution(this.mean, this.cv);
		this.useRootEdge = useRootEdge;
		this.logLikelihoods = null;
	}

	// Update the log likelihoods
	public void update(Tree tree) {
		this.length = tree.bl;
		this.time = new double[tree.nnodes];
		this.logLikelihoods = new LogDouble[tree.nnodes];
		for (Node node : tree.nodes) {
			// Just a special case: we do not consider the rate model over the stem arc if
			// desired (when doUseRootEdge = false).
			boolean doUseRootEdge = (this.useRootEdge || !node.isRoot());
			if (doUseRootEdge) {
				time[node.id] = node.isRoot() ? tree.peakTime - tree.vt[node.id]
						: tree.vt[node.parent.id] - tree.vt[node.id];
				if (time[node.id] != 0.0) {
					double rate = this.length[node.id] / this.time[node.id];
					logLikelihoods[node.id] = edgeRatePD.getPDFAsProbability(rate);
				}
			}

		}
	}

	// Get complete model likelihood
	public LogDouble getLogLikelihood() {
		LogDouble modelLikelihood = new LogDouble(1.0);
		if (logLikelihoods != null) {
			for (int i = 0; i < logLikelihoods.length; i++) {
				if (logLikelihoods[i] != null)
					modelLikelihood.mult(logLikelihoods[i]);
			}
		}
		return modelLikelihood;
	}

	// Getter method for Gamma distribution
	public double getMean() {
		return this.mean;
	}

	public double getCV() {
		return this.cv;
	}

	// Show the gamma parameters
	public void ShowModelParameters() {
		System.out.println("Gamma Dist. mean =" + this.edgeRatePD.getMean());
		System.out.println("Gamma Dist. CV   =" + this.edgeRatePD.getCV());
		System.out.println("Gamma k parameter=" + this.edgeRatePD.k);
		System.out.println("Gamma theta parameter=" + this.edgeRatePD.theta);
	}

}

/*
 * // First update the log likelihoods for each node then return model
 * likelihood public LogDouble getUpdatedLogLikelihood() {
 * 
 * this.completeUpdate();
 * 
 * LogDouble modelLikelihood = new LogDouble(1.0); for(Node node: tree.nodes) {
 * if( !node.isRoot() ) modelLikelihood.mult(logLikelihoods[node.id]); } return
 * modelLikelihood; }
 * 
 * //This method will serve local likelihood improvement along a node public
 * LogDouble getPartialLogLikelihood(Node node) {
 * 
 * double time,length; LogDouble patialLogLikelihood = new LogDouble(1.0);
 * 
 * if(node.isRoot()) { for(int i=0 ; i<node.nchildren ; i++) { Node child =
 * node.children.get(i); time = tree.vt[child.parent.id] - tree.vt[child.id];
 * length = tree.bl[child.id]; patialLogLikelihood.mult(
 * edgeRatePD.getPDFAsProbability(length/time) ) ; } return patialLogLikelihood;
 * }
 * 
 * 
 * time = tree.vt[node.parent.id] - tree.vt[node.id]; length = tree.bl[node.id];
 * patialLogLikelihood.mult( edgeRatePD.getPDFAsProbability(length/time) );
 * for(int i=0 ; i<node.nchildren ; i++) { Node child = node.children.get(i);
 * time = tree.vt[child.parent.id] - tree.vt[child.id]; length =
 * tree.bl[child.id]; patialLogLikelihood.mult(
 * edgeRatePD.getPDFAsProbability(length/time) ) ; }
 * 
 * 
 * return patialLogLikelihood; }
 * 
 * 
 * public static LogDouble getLogLikelihood1(double [] edgeRates , Tree tree,
 * Log log) throws IOException {
 * 
 * double mean = edgeRates[0]; double cv = edgeRates[1];
 * 
 * double [] lengths = tree.bl; double [] arcTimes = tree.at;
 * 
 * LogDouble modelLikelihood = new LogDouble(1.0); // Model likelihood LogDouble
 * [] logLikelihoods = new LogDouble [tree.nnodes]; // Node wise likelihood
 * GammaDistribution edgeRatePD = new GammaDistribution(mean,cv); // Probability
 * distribution StringBuilder logStr = new StringBuilder(); // String buffer for
 * logging
 * 
 * for(Node node: tree.nodes) { if( !node.isRoot() ) {
 * 
 * logLikelihoods[node.id] = edgeRatePD.getPDFAsProbability( lengths[node.id] /
 * arcTimes[node.id] ) ; if(Double.isNaN(logLikelihoods[node.id].getLogValue()))
 * { logStr.append("\nRate Model Error: NaN issue with the");
 * logStr.append(tree.name + "\t" + "Node id = " + node.id + "\tlen =" +
 * lengths[node.id] + "\ttime =" + arcTimes[node.id]); continue; }
 * modelLikelihood.mult(logLikelihoods[node.id]); } } if(log!=null)
 * log.write(logStr.toString());
 * 
 * return modelLikelihood; }
 */
