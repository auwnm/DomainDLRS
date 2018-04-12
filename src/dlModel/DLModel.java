package dlModel;

import java.io.IOException;
import java.util.ArrayList;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.analysis.solvers.UnivariateSolver;

import common.Log;
import common.LogDouble;
import common.Mathematics;
import common.PRNG;
import phylogeny.Node;
import phylogeny.Reconciliation;
import phylogeny.Tree;

//Adapted from Matt Rasmussen 2011 Gene tree topology prior of the SPIMAP model
public class DLModel {

	private Log log;
	private PRNG prng;
	private Reconciliation R;
	private double[] BDRates;
	private double[] extinctionProb;
	private boolean useRootEdge;

	public DLModel(Reconciliation R, double[] BDRates, double[] extinctionProb, boolean useRootEdge, PRNG prng,
			Log log) {
		this.R = R;
		this.BDRates = BDRates;
		this.extinctionProb = extinctionProb;
		this.useRootEdge = useRootEdge;
		this.prng = prng;
		this.log = log;
	}

	public void updateDLModel(double[] extinctionProb) {
		this.extinctionProb = extinctionProb;
	}

	// Returns the number of "Equally Likely Labeled" Histories for 'ngenes'
	// survivors
	private double ellHistories(int ngenes) {
		double n = 1;
		for (int i = 2; i <= ngenes; i++) {
			n *= i * (i - 1) / 2;
		}
		return n;
	}

	// Returns the number of labeled histories exist for the given tree topology
	// uses the subtree starting at root and going until leaves.
	// NOTE: assumes binary tree
	private double numSubtopologyHistories(Tree tree, Node root, ArrayList<Node> leaves) {
		double n = 1.0;

		// get nodes in post order
		ArrayList<Node> queue = new ArrayList<Node>(tree.nnodes);

		// helper array for ensuring postorder traversal
		int[] visited = new int[tree.nnodes];
		for (int i = 0; i < visited.length; i++)
			visited[i] = 0;

		// count number of descendant internal nodes
		int[] ninternals = new int[tree.nnodes];

		// process leaves
		for (int i = 0; i < leaves.size(); i++) {
			Node node = leaves.get(i);
			ninternals[node.id] = 0;

			// queue parent
			visited[node.parent.id]++;
			queue.add(node.parent);
		}

		// go up tree until root
		for (int i = 0; i < queue.size(); i++) {
			Node node = queue.get(i);

			// do not process a node until both children are processed
			if (visited[node.id] != 2)
				continue;

			// count internal children
			int right = ninternals[node.children.get(1).id];
			int left = ninternals[node.children.get(0).id];
			ninternals[node.id] = 1 + right + left;
			n *= Mathematics.fchoose(right + left, right);

			if (node == root)
				return n;

			visited[node.id]++;
			visited[node.parent.id]++;
			queue.add(node.parent);
		}

		// Note: this will occur for subtree like
		// root
		// |
		// +
		// / \
		// / \
		// A B
		//
		return n;
	}

	// Get the leaves of the subtree below gnode that reconciles to snode
	private void getSubTreeLeaves(Node gnode, Node snode, ArrayList<Node> leaves) {
		if (leaves == null)
			leaves = new ArrayList<Node>();

		for (int i = 0; i < gnode.nchildren; i++) {
			Node child = gnode.children.get(i);
			// Only consider nodes that reconcile to snode
			if (this.R.nodeMapPrime[child.id] == snode.id) {
				if (this.R.eventMapPrime[child.id] == Reconciliation.EVENT_SPEC
						|| this.R.eventMapPrime[child.id] == Reconciliation.EVENT_LEAF)
					leaves.add(child);
				else // duplication case
					getSubTreeLeaves(child, snode, leaves);
			}
		}
	}

	// Get the duplications of the subtree below gnode that reconciles to snode
	private void getSubTreeDups(Node gnode, Node snode, ArrayList<Node> duplications) {
		if (duplications == null)
			duplications = new ArrayList<Node>();

		for (int i = 0; i < gnode.nchildren; i++) {
			Node child = gnode.children.get(i);
			// Only consider internal nodes that reconcile to this snode
			if (this.R.nodeMapPrime[child.id] == snode.id) {
				if (this.R.eventMapPrime[child.id] == Reconciliation.EVENT_DUPL) {
					duplications.add(child);
					getSubTreeDups(child, snode, duplications);
				}
			}
		}
	}

	// Probability of Subtree
	// gnode
	// |
	// +
	// / \
	// + +
	// / \ / \
	// Leaves
	private LogDouble getSubTreeProb(Node gnode, Node snode, double t, double extinctionProb) {

		double birthRate = this.BDRates[0];
		double deathRate = this.BDRates[1];

		// Initializing sub-tree variables
		double phist;
		double subTreeProb = 0.0;

		// Getting leaves of the subtree reconcile to this snode
		ArrayList<Node> leaves = new ArrayList<Node>();
		getSubTreeLeaves(gnode, snode, leaves);
		int s = leaves.size();

		// Computing P0, P1 and denominator term a
		double r = birthRate / deathRate;
		double p0 = BirthDeathProbs.P0(birthRate, deathRate, t);
		double p1 = BirthDeathProbs.P1(birthRate, deathRate, t);
		double a = 1.0 - r * p0 * extinctionProb;

		// Evaluating Cases
		if (s > 1) {
			phist = numSubtopologyHistories(this.R.guestPrime, gnode, leaves) / ellHistories(s);
			subTreeProb = phist * Math.pow(r * p0, s - 1) * p1 / Math.pow(a, s + 1);
		} else if (s == 1) {
			phist = numSubtopologyHistories(this.R.guestPrime, gnode, leaves) / ellHistories(s);
			subTreeProb = (phist * p1) / (a * a);
		} else if (s == 0) {
			subTreeProb = p0 + (p1 * extinctionProb) / a;
		}

		LogDouble probability = new LogDouble(subTreeProb);
		return probability;
	}

	private boolean UltrametrizeSubTree(Node gnode, Node snode, double t, double dy) throws IOException {

		// Getting leaves of the subtree reconcile to this snode
		ArrayList<Node> leaves = new ArrayList<Node>();
		getSubTreeLeaves(gnode, snode, leaves);

		// Setting vertex times of leaves according to species vertex
		for (Node node : leaves)
			this.R.guestPrime.vt[node.id] = R.host.vt[snode.id];

		// Getting duplications of the subtree reconcile to this snode
		ArrayList<Node> duplications = new ArrayList<Node>();
		getSubTreeDups(gnode, snode, duplications);
		if (duplications.size() == 0)
			return true;

		// Sample time for each duplication (internal node)
		boolean correctlySampled = false;
		correctlySampled = sampleTime(duplications.get(0), snode, dy, t);

		return correctlySampled;
	}

	// Recursive call to sample the times for internal nodes of this subtree
	private boolean sampleTime(Node node, Node snode, double deathProbAtY, double t) throws IOException {

		double sampledTime = 0;
		double birthRate = this.BDRates[0];
		double deathRate = this.BDRates[1];
		double[] vertexTimes = this.R.guestPrime.vt;

		// Initializing sub-tree variables
		ArrayList<Node> subleaves = new ArrayList<Node>();
		ArrayList<Node> subduplications = new ArrayList<Node>();

		getSubTreeLeaves(node, snode, subleaves);
		getSubTreeDups(node, snode, subduplications);

		int s = subleaves.size();
		double phist = numSubtopologyHistories(this.R.guestPrime, node, subleaves) / ellHistories(s);
		double rValue = prng.nextDouble();
		double cdfFunctionValue = Math.log(rValue);

		// Setting Numerical analysis solver parameters
		double minValue = 0;
		double maxValue = t;
		int maxIterations = 10000;
		double relativeAccuracy = 1.0e-12;
		double absoluteAccuracy = 1.0e-8;
		TimeSampler ts = new TimeSampler(birthRate, deathRate, t, deathProbAtY, phist, s, cdfFunctionValue);

		if (t == 0.0) {
			log.write("\nTimeSampler issue: Time span for sampling becomes zero.");
			return false;
		}
		if (ts.SubTreeProbIsZero()) {
			log.write("\nTimeSampler issue: Subtree probability approaches to zero as time = " + t);
			return false;
		}
		double fmin = ts.value(minValue);
		double fmax = ts.value(maxValue);
		if ((fmin * fmax) >= 0.0) {
			log.write("\nTimeSampler issue: Root is not bracketd within specified interval.");
			return false;
		}

		// Root finding
		double root = 0;
		UnivariateSolver nonBracketing = new BrentSolver(relativeAccuracy, absoluteAccuracy);
		try {
			root = nonBracketing.solve(maxIterations, ts, minValue, maxValue);
		} catch (RuntimeException re) {
			log.newLine();
			log.write("TimeSampler issue: Numerical solver failed to find the root.");
			return false;
		}
		sampledTime = t - root;
		vertexTimes[node.id] = R.host.vt[snode.id] + sampledTime;

		// Recursion for both childrens (assume binary tree)
		boolean correctlySampled = true;
		for (int i = 0; i < node.children.size(); i++) {
			Node child = node.children.get(i);
			if (this.R.nodeMapPrime[child.id] == snode.id
					&& this.R.eventMapPrime[child.id] == Reconciliation.EVENT_DUPL)
				correctlySampled &= sampleTime(child, snode, deathProbAtY, sampledTime);
		}

		return correctlySampled;
	}

	public String showSubTreeStr(Node node, Node snode) {

		ArrayList<Node> leaves = new ArrayList<Node>();
		ArrayList<Node> duplications = new ArrayList<Node>();
		getSubTreeLeaves(node, snode, leaves);
		getSubTreeDups(node, snode, duplications);

		StringBuilder sb = new StringBuilder();
		sb.append("(" + snode.id + " -< ");

		if (duplications.size() != 0) {
			sb.append("(" + snode.id + " -< ");
			for (int i = 0; i < duplications.size() - 1; i++)
				sb.append(duplications.get(i).id + ",");
			sb.append(duplications.get(duplications.size() - 1).id + " : ");

		} else
			sb.append("(" + snode.id + " : ");

		for (int i = 0; i < leaves.size() - 1; i++)
			sb.append(leaves.get(i).id + ",");
		sb.append(leaves.get(leaves.size() - 1).id + ") ");

		return sb.toString();
	}

	public boolean sampleTimes() throws IOException {

		Tree guest = this.R.guestPrime;
		Tree host = this.R.host;
		int[] nodeMap = this.R.nodeMapPrime;
		int[] eventMap = this.R.eventMapPrime;

		boolean correctlySampled = true;

		// Initialize the time arrays of guest tree
		guest.vt = new double[guest.nnodes];

		// Handling preroot duplications
		if (eventMap[guest.root.id] == Reconciliation.EVENT_DUPL) {

			// Adding pseudo root
			Node pseudoRoot = new Node();
			guest.addNode(pseudoRoot);
			pseudoRoot.addChild(guest.root);
			guest.root.parent = pseudoRoot;
			guest.root = pseudoRoot;

			Node snode = host.root;
			double t = host.getPeakTime() - host.vt[snode.id];
			double dc = Math.exp(extinctionProb[snode.id]);

			correctlySampled &= UltrametrizeSubTree(pseudoRoot, snode, t, dc);

			// Removing pseudo root
			guest.root = pseudoRoot.children.get(0);
			guest.root.parent = null;
			guest.nodes.remove(guest.nnodes - 1);
			guest.nnodes--;
		}

		// Setting time of root node
		if (eventMap[guest.root.id] == Reconciliation.EVENT_SPEC) {
			assert (R.nodeMap[guest.root.id] == R.host.root.id);
			guest.vt[guest.root.id] = host.vt[host.root.id];
		}

		// Loop through speciation nodes in tree
		for (int i = 0; i < guest.nnodes; i++) {
			Node node = guest.nodes.get(i);
			if (eventMap[node.id] == Reconciliation.EVENT_SPEC) {
				// loop through nodes u \in child(R(v))
				Node snode = host.nodes.get(nodeMap[node.id]);
				for (int j = 0; j < snode.nchildren; j++) {

					// Getting the subtree probability
					Node schild = snode.children.get(j);
					double t = host.vt[schild.parent.id] - host.vt[schild.id];
					double dc = Math.exp(extinctionProb[schild.id]);

					correctlySampled &= UltrametrizeSubTree(node, schild, t, dc);

				}
			} else if (eventMap[node.id] == Reconciliation.EVENT_DUPL) {
				// modelLikelihood.mult(Math.log(2));
			}

		}

		return correctlySampled;

	}

	// TODO: does not handle branches above the species tree root yet
	// NOTE: assumes binary species tree
	public boolean sampleRealisation() throws IOException {

		// Sample times for guest tree in Reconciliation
		boolean correctlySampled = sampleTimes();

		// This method will copy vertex times from implied vertices tree
		// i.e. for first n vertices also set peak time for root node
		if (R.guest.vt == null)
			R.guest.vt = new double[R.guest.nnodes];
		for (int idx = 0; idx < R.guest.nnodes; idx++) {
			Node node = R.guest.nodes.get(idx);
			R.guest.vt[node.id] = R.guestPrime.vt[node.id];
		}
		R.guest.peakTime = R.host.peakTime;

		if (!correctlySampled)
			return false;
		else
			return isValidRealization(this.R);
	}

	// This method will validate the realization
	public boolean isValidRealization(Reconciliation R) throws IOException {

		boolean isValid = true;
		StringBuilder sb = new StringBuilder();

		// Check vertex time for each node
		for (Node node : R.guest.nodes) {

			double vertexTime = R.guest.vt[node.id];
			int hostIdx = R.nodeMap[node.id];

			if (vertexTime < 0.0) {
				sb.append("\nInvalid realization of " + R.guest.name + " :\t" + "vertex time of node." + node.id
						+ " less than zero.");
				isValid = false;
			}
			if (vertexTime != 0.0 && R.eventMap[node.id] == Reconciliation.EVENT_LEAF) {
				sb.append("\nInvalid realization of " + R.guest.name + " :\t" + "vertex time of leaf." + node.id
						+ " is not equal to zero.");
				isValid = false;
			} else if (vertexTime != R.host.vt[hostIdx] && R.eventMap[node.id] == Reconciliation.EVENT_SPEC) {
				sb.append("\nInvalid realization of " + R.guest.name + " :\t" + "vertex time of node." + node.id
						+ " is not equal to host speciation node time.");
				isValid = false;
			} else if (R.eventMap[node.id] == Reconciliation.EVENT_DUPL) {
				double upperTime;
				double lowerTime = R.host.vt[hostIdx];
				if (R.host.nodes.get(hostIdx).isRoot()) {
					upperTime = R.host.getPeakTime();
				} else {
					int parentIdx = R.host.nodes.get(hostIdx).parent.id;
					upperTime = R.host.vt[parentIdx];
				}
				if (!(R.guest.vt[node.id] > lowerTime && R.guest.vt[node.id] < upperTime)) {
					isValid = false;
					sb.append("\nInvalid realization of " + R.guest.name + ":\t" + "vertex time of duplication node."
							+ node.id + " time = " + vertexTime + " does not lie within " + "[ " + upperTime + " , "
							+ lowerTime + " ]");
				}
			}
		}

		// Writing to log
		if (sb.length() != 0) {
			log.newLine();
			log.write(sb.toString());
		}

		return isValid;
	}

	// TODO: does not handle branches above the species tree root yet
	// NOTE: assumes binary species tree
	public LogDouble birthDeathTreePrior() throws IOException {

		double birthRate = this.BDRates[0];
		double deathRate = this.BDRates[1];

		Tree guest = this.R.guestPrime;
		Tree host = this.R.host;
		int[] nodeMap = this.R.nodeMapPrime;
		int[] eventMap = this.R.eventMapPrime;

		LogDouble subTreeProb;
		LogDouble modelLikelihood = new LogDouble(1.0);

		// Handling preroot duplications
		if (eventMap[guest.root.id] == Reconciliation.EVENT_DUPL) {

			// Adding pseudo root
			Node pseudoRoot = new Node();
			guest.addNode(pseudoRoot);
			pseudoRoot.addChild(guest.root);
			guest.root.parent = pseudoRoot;
			guest.root = pseudoRoot;

			Node snode = host.root;
			double t = host.getPeakTime() - host.vt[snode.id];
			double dc = Math.exp(extinctionProb[snode.id]);

			subTreeProb = getSubTreeProb(pseudoRoot, host.root, t, dc);
			modelLikelihood.mult(subTreeProb);

			// Removing pseudo root
			guest.root = pseudoRoot.children.get(0);
			guest.root.parent = null;
			guest.nodes.remove(guest.nnodes - 1);
			guest.nnodes--;
		}

		// Computing subTree probability of stem edge i.e. P11 for single lineage
		if (eventMap[guest.root.id] == Reconciliation.EVENT_SPEC && this.useRootEdge) {
			assert (R.nodeMap[guest.root.id] == R.host.root.id);

			Node snode = host.root;
			double t = host.getPeakTime() - host.vt[snode.id];
			double dc = Math.exp(extinctionProb[snode.id]);
			double p0 = BirthDeathProbs.P0(birthRate, deathRate, t);
			double p1 = BirthDeathProbs.P1(birthRate, deathRate, t);
			double a = 1.0 - (birthRate / deathRate) * p0 * dc;
			subTreeProb = new LogDouble((p1 / a * a));
			modelLikelihood.mult(subTreeProb);
		}

		// Loop through speciation nodes in tree
		for (int i = 0; i < guest.nnodes; i++) {
			Node node = guest.nodes.get(i);
			if (eventMap[node.id] == Reconciliation.EVENT_SPEC) {
				// loop through nodes u \in child(R(v))
				Node snode = host.nodes.get(nodeMap[node.id]);
				for (int j = 0; j < snode.nchildren; j++) {

					// Getting the subtree probability
					Node schild = snode.children.get(j);
					double t = host.vt[schild.parent.id] - host.vt[schild.id];
					double dc = Math.exp(extinctionProb[schild.id]);
					subTreeProb = getSubTreeProb(node, schild, t, dc);
					modelLikelihood.mult(subTreeProb);

				}
			} else if (eventMap[node.id] == Reconciliation.EVENT_DUPL) {
				// modelLikelihood.mult(Math.log(2));
			}

		}
		return modelLikelihood;
	}

}

// double t = host.at[schild.id];
// double t = host.at[schild.id];
// public static boolean firstIteration = true;

/*
 * if(snode.isRoot() && this.R.guest.name.equalsIgnoreCase("domainTree2") ) {
 * System.out.println("--------------------"); for(int i=0 ;
 * i<duplications.size() ; i++) { Node dup = duplications.get(i);
 * System.out.println( "Node " + dup.id + " time: \t " +
 * this.R.guestPrime.vt[dup.id]); } System.out.println("Lower time = " +
 * this.R.host.vt[snode.id]); System.out.println("--------------------"); }
 * 
 * if(snode.isRoot() && this.R.guest.name.equalsIgnoreCase("domainTree2") &&
 * firstIteration ) { System.out.println("Number of duplications = " +
 * duplications.size()); System.out.println("Number of Leaves = " +
 * leaves.size()); System.out.println("Birth rate = " + this.BDRates[0]);
 * System.out.println("Death rate = " + this.BDRates[1]);
 * System.out.println("Pe at Y = " + dy);
 * System.out.println("Time Iterval of gene stem = " + "[" +
 * this.R.host.vt[snode.id] + " , " + (this.R.host.vt[snode.id] + t) +"]");
 * System.out.println("Leaves :" ); for(Node leaf: leaves)
 * System.out.print(leaf.id + ", "); System.out.println();
 * System.out.println("----------------------------------------------------");
 * isSampled = sampleTime(duplications.get(0),snode,dy,t,prng,true);
 * System.out.println("----------------------------------------------------");
 * firstIteration = false; } else
 */

/*
 * public static double MIN_BRANCH_TIME = 1e-3; // 0.000001 public static int
 * MAX_REAL_TRIES = 10; // 100 int tries = 0; boolean isValid = false; boolean
 * isSampled = false; do { isSampled =
 * sampleTime(duplications.get(0),snode,dy,t,prng,false);
 * 
 * for(int i=0 ; i<duplications.size() ; i++) {
 * 
 * Node dup = duplications.get(i);
 * 
 * // time to parent double parentTime = (i==0) ? this.R.host.vt[snode.id] + t :
 * this.R.guestPrime.vt[dup.parent.id] ; double t0 = parentTime -
 * this.R.guestPrime.vt[dup.id] ;
 * 
 * // time to children double t1 = this.R.guestPrime.vt[dup.id] -
 * this.R.guestPrime.vt[dup.children.get(0).id]; double t2 =
 * this.R.guestPrime.vt[dup.id] - this.R.guestPrime.vt[dup.children.get(1).id];
 * 
 * if( t0 < MIN_BRANCH_TIME || t1 < MIN_BRANCH_TIME || t2 < MIN_BRANCH_TIME )
 * isValid = false;
 * 
 * } if(tries ++ > MAX_REAL_TRIES) break;
 * 
 * }while( !isValid && !isSampled );
 */

// For Debuging
// if(node.id == 10 && this.R.guest.name.equalsIgnoreCase("geneTree"))
// System.out.println(vertexTimes[node.id] + "\t" + rValue );
/*
 * if(snode.isRoot() && this.R.guest.name.equalsIgnoreCase("domainTree2") &&
 * debug) { System.out.println( "duplication." + node.id +" parent of ("+
 * node.children.get(0).id +" , "+ node.children.get(1).id + ")" +
 * " sampled time: \t " + this.R.guestPrime.vt[node.id] if(node.id == 128 ) {
 * TimeSampler ts1 = new TimeSampler(birthRate, deathRate, t, deathProbAtY,
 * phist, s, cdfFunctionValue ); double delta = t/1000; double to = 0.0; do {
 * System.out.println( to + "\t" + ts1.value1(to) ); to += delta;
 * 
 * }while(to<=t); } }
 */

////////////////////////////////////////////////
/*
 * // This method will check the arc times if these are very small ... boolean
 * tooSmallArcTimes(Reconciliation R) { double arcTime; boolean isTooSmall =
 * false;
 * 
 * for(Node node: R.guest.nodes) { if( node.isRoot() ) arcTime =
 * R.guest.getPeakTime() - R.guest.vt[node.id]; else arcTime =
 * R.guest.vt[node.parent.id] - R.guest.vt[node.id];
 * 
 * if( arcTime < MIN_BRANCH_TIME ) isTooSmall = true; } return isTooSmall; }
 * 
 * if(snode.isRoot() && this.R.guest.name.equalsIgnoreCase("domainTree2") ) {
 * System.out.println("--------------------"); for(int i=0 ;
 * i<duplications.size() ; i++) { Node dup = duplications.get(i);
 * System.out.println( "Node " + dup.id + " time: \t " +
 * this.R.guestPrime.vt[dup.id]); } System.out.println("Lower time = " +
 * this.R.host.vt[snode.id]); System.out.println("--------------------"); }
 * 
 * 
 * if(snode.isRoot() && this.R.guest.name.equalsIgnoreCase("geneTree") ) {
 * System.out.println("pre root gene tree duplication .... ");
 * 
 * }
 */

// Convenience functions
// Update time from implied species nodes to the gene tree
// NOTE: Assumes binary species tree & tree with implied nodes
/*
 * public static LogDouble getBirthDeathTreePriorFull(Reconciliation R, double
 * [] BDRates,double [] extinctionProb,Log log) throws IOException { LogDouble
 * loglikelihood = birthDeathTreePrior(R,BDRates ,extinctionProb,true,true,log);
 * updateSampledTime(R,log); return loglikelihood; }
 * 
 * 
 * public static LogDouble getBirthDeathTreePrior(Reconciliation R, double []
 * BDRates , double [] extinctionProb,Log log) throws IOException { return
 * birthDeathTreePrior(R,BDRates ,extinctionProb,true,false,log); }
 * 
 * 
 * 
 * public static void sampleRealisation(Reconciliation R, double [] BDRates ,
 * double [] extinctionProb,Log log) throws IOException {
 * birthDeathTreePrior(R,BDRates ,extinctionProb,false,true,log);
 * updateSampledTime(R,log); return; }
 * 
 * 
 * // Setting arc time from vertex time // R.guest.setArcTimesFromVertexTimes();
 * // Setting arc time for root node
 * 
 * int rid = R.host.root.id; double peakTime = R.host.vt[rid] + R.host.at[rid];
 * R.guest.at[R.guest.root.id] = peakTime - R.guestPrime.vt[R.guest.root.id];
 * 
 * // Checking realization ultrametric property if(!R.guest.isUltrametric())
 * log.write("\nDL Model Error: Non ultrametric issue with the " +
 * R.guest.name);
 * 
 */
