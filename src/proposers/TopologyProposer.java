package proposers;

import java.util.ArrayList;

import Main.Parameters;
import phylogeny.Node;
import phylogeny.Reconciliation;
import phylogeny.Tree;

import common.LogDouble;
import common.Mathematics;
import common.PRNG;
import common.Pair;
import common.RealInterval;

public class TopologyProposer implements Proposer {

	public String name;
	public Parameters pi;
	public int level;
	public int domIndex;

	public PRNG prng;
	public double[] moveWeights;
	public String lastOperationType;
	public RealInterval interval;
	public ProposerStatistics stats;

	public Tree treeCache; // Tree parameter cache.
	public Pair<Node, Node> nodePair; // Pair of nodes to specify the NNI,SPR,REROOT location
	public Node sprNodec; // For remembering the sibling of nodea in case of SPR

	private LogDouble forwardDensity; // Forward probability density.
	private LogDouble backwardDensity; // Backward probability density.
	private boolean isEnabled; // On/off switch. By default it will be true.

	public TopologyProposer(Parameters pi, int level, int domIndex, double[] moveWeights, PRNG prng) {
		this.pi = pi;
		this.level = level;
		this.domIndex = domIndex;

		if (this.level == Parameters.GENE_LEVEL)
			this.name = "geneTreeTopology";
		else
			this.name = "domTreeTopologyProp_" + (this.domIndex + 1);

		this.treeCache = null;
		this.nodePair = null;
		this.sprNodec = null;
		this.prng = prng;
		this.isEnabled = true;
		this.moveWeights = moveWeights;
		this.stats = new ProposerStatistics(this.name);
	}

	/*
	 * public static void main(String[] args) {
	 * 
	 * int rv;
	 * 
	 * rv= Mathematics.irand(5); System.out.println(rv);
	 * 
	 * rv = Mathematics.irand(5); System.out.println(rv);
	 * 
	 * }
	 */

	// This method will perturb branch lengths w.r.t. its tree structure
	// nextInt = Returns a pseudorandom, uniformly distributed int value between 0
	// (inclusive) and the specified value (exclusive)
	@Override
	public boolean cacheAndPerturb() {

		Tree tree;
		Reconciliation recon;

		// First specify the tree paramter
		if (this.level == Parameters.GENE_LEVEL) {
			tree = pi.geneTree;
			recon = new Reconciliation(pi.geneTree, pi.speciesTree, pi.gene2species);
		} else {
			tree = pi.domainTree[this.domIndex];
			recon = new Reconciliation(pi.domainTree[this.domIndex], pi.geneTree, pi.domain2genes[this.domIndex]);
		}

		// Second maintain the current state of tree parameter
		tree.cacheLengths();
		this.treeCache = tree.copy();

		// Perturb: First determine move to make.
		double w = this.prng.nextDouble()
				* (this.moveWeights[0] + this.moveWeights[1] + this.moveWeights[2] + this.moveWeights[3]);
		if (w < this.moveWeights[0]) {
			this.doNNI(tree);
			this.lastOperationType = "NNI";
		} else if (w < this.moveWeights[0] + this.moveWeights[1]) {
			this.doRNNI(recon, this.prng);
			this.lastOperationType = "RNNI";
		} else if (w < this.moveWeights[0] + this.moveWeights[1] + this.moveWeights[2]) {
			this.doSPR(tree);
			this.lastOperationType = "SPR";
		} else if (w < this.moveWeights[0] + this.moveWeights[1] + this.moveWeights[2] + this.moveWeights[3]) {
			this.doREROOTING(tree);
			this.lastOperationType = "REROOTING";
		}

		/*
		 * if(this.level == Parameters.GENE_LEVEL)
		 * System.out.println(this.lastOperationType);
		 */

		// Set forward and backward densities (we consider forward-backward
		// probabilities as equal)
		this.forwardDensity = new LogDouble(1.0);
		this.backwardDensity = new LogDouble(1.0);

		tree.updateTree();

		return true;
	}

	// Clears the cached (previous) state of perturbed parameters when a proposed
	// state has been accepted.
	@Override
	public void clearCache() {
		if (this.stats != null) {
			this.stats.increment(true, this.lastOperationType);
		}

		Tree currentTree;
		// First specify the tree paramter
		if (this.level == Parameters.GENE_LEVEL)
			currentTree = pi.geneTree;
		else
			currentTree = pi.domainTree[this.domIndex];

		currentTree.clearLengthsCache();

		this.nodePair = null;
		this.sprNodec = null;
		this.treeCache = null;
	}

	// Restores the cached (previous) state of perturbed parameters when a proposed
	// state has been rejected.
	@Override
	public void restoreCache() {
		if (this.stats != null) {
			this.stats.increment(false, this.lastOperationType);
		}

		Tree currentTree;

		// First specify the tree paramter
		if (this.level == Parameters.GENE_LEVEL)
			currentTree = pi.geneTree;
		else
			currentTree = pi.domainTree[this.domIndex];

		// Recovering length part
		if (currentTree.lengthCache != null) {
			currentTree.bl = currentTree.lengthCache;
			this.treeCache.bl = currentTree.lengthCache;
			currentTree.clearLengthsCache();
		} else {
			this.treeCache.bl = currentTree.bl;
		}

		// Replacing with new tree
		if (this.level == Parameters.GENE_LEVEL) {
			pi.geneTree = this.treeCache;
			pi.geneTree.updateTree();
		} else {
			pi.domainTree[this.domIndex] = this.treeCache;
			pi.domainTree[this.domIndex].updateTree();
		}

		this.sprNodec = null;
		this.nodePair = null;
		this.treeCache = null;

	}

	// Propose n Perform NNI & SPR
	private boolean doNNI(Tree tree) {
		this.nodePair = proposeRandomNni(tree);
		performNni(tree, this.nodePair);
		return true;
	}

	// Propose n Perform reconciliation based NNI
	private boolean doRNNI(Reconciliation recon, PRNG prng) {
		this.nodePair = null;
		this.nodePair = performRNNI(recon, prng);

		if (this.nodePair == null)
			return false;

		return true;
	}

	private boolean doSPR(Tree tree) {
		this.nodePair = proposeRandomSpr(tree);
		Node pnode = nodePair.first.parent;
		this.sprNodec = (pnode.children.get(0) == this.nodePair.first) ? pnode.children.get(1) : pnode.children.get(0);

		if (!validSpr(tree, this.nodePair.first, this.nodePair.second))
			return false;

		performSpr(tree, this.nodePair.first, this.nodePair.second);
		return true;
	}

	private boolean doREROOTING(Tree tree) {

		// First save old root position
		Node oldRoot1 = tree.root.children.get(0);
		Node oldRoot2 = tree.root.children.get(1);
		this.nodePair = new Pair<Node, Node>(oldRoot1, oldRoot2);

		Pair<Node, Node> edge = proposeRandomRerooting(tree);
		performRerooting(tree, edge);

		return true;
	}

	// Returns the "forward" probability density Q(x';x) for obtaining the new value
	// x' given the old value x.
	@Override
	public LogDouble getForwardDensity() {
		return this.forwardDensity;
	}

	// Returns the "backward" probability density Q(x;x') for obtaining the old
	// value x given the new value x'.
	@Override
	public LogDouble getBackwardDensity() {
		return this.backwardDensity;
	}

	// The ratio between the "backward" and "forward" proposal densities
	// respectively.
	// It will returns the ratio Q(x;x')/Q(x';x) for the old state x and the new
	// state x',
	@Override
	public LogDouble getDensityRatio() {
		return this.backwardDensity.divToNew(this.forwardDensity);
	}

	// This method will return the name of acting proposer
	@Override
	public String getProposerName() {
		return this.name;
	}

	// This method will return the object of statistics associated with this
	// proposer
	@Override
	public ProposerStatistics getProposerStatistics() {
		return this.stats;
	}

	@Override
	public boolean isEnabled() {
		return this.isEnabled;
	}

	@Override
	public void setEnabled(boolean isActive) {
		this.isEnabled = isActive;
	}

	@Override
	public boolean hasValidProposal() {
		return (!this.backwardDensity.isZero() && !this.forwardDensity.isZero());
	}

	// This will return the index associated with perturbed rates
	public String getOperationType() {
		return this.lastOperationType;
	}

	public void setMoveWeights(double[] moveWeights) {
		this.moveWeights = moveWeights;
	}

	// =============================================================================
	// Nearest Neighbor Interchange Topology Proposal
	// Derived from Rasmussen et. al. 2011
	/*
	 * 
	 * Proposes a new tree using Nearest Neighbor Interchange
	 * 
	 * Branch for NNI is specified by giving its two incident nodes (node1 and
	 * node2). Change specifies which subtree of node1 will be swapped with the
	 * uncle. See figure below.
	 * 
	 * node2 / \ nodeb node1 / \ nodea *
	 * 
	 */

	public static void performNni(Tree tree, Pair<Node, Node> nniPair) {

		Node nodea = nniPair.first;
		Node nodeb = nniPair.second;

		Node node1 = nodea.parent;
		Node node2 = nodeb.parent;

		// assert that node1 and node2 are incident to the same branch
		assert (node1.parent == node2 || node2.parent == node1);

		// find child indexes
		int a = (node1.children.get(0) == nodea) ? 0 : 1;
		assert (node1.children.get(a) == nodea);

		int b = (node2.children.get(0) == nodeb) ? 0 : 1;
		assert (node2.children.get(b) == nodeb);

		// swap parent pointers
		nodea.parent = node2;
		nodeb.parent = node1;

		// swap child pointers
		node2.children.set(b, nodea);
		node1.children.set(a, nodeb);
	}

	public static Pair<Node, Node> proposeRandomNni(Tree tree) {
		Node a, b;

		// find edges for NNI
		int choice;
		do {
			choice = Mathematics.irand(tree.nnodes);
		} while (tree.nodes.get(choice).isLeaf() || tree.nodes.get(choice).isRoot());

		Node node1 = tree.nodes.get(choice);
		Node node2 = tree.nodes.get(choice).parent;
		a = node1.children.get(Mathematics.irand(2));
		b = (node2.children.get(0).id == node1.id) ? node2.children.get(1) : node2.children.get(0);

		assert (a.parent.parent == b.parent);

		Pair<Node, Node> nniPair = new Pair<Node, Node>(a, b);
		return nniPair;
	}

	public static Pair<Node, Node> performRNNI(Reconciliation recon, PRNG prng) {

		int smallest_nnodes = 20;
		Tree tree = recon.guest;

		if (tree.nnodes < smallest_nnodes) {
			Pair<Node, Node> nniPair;
			nniPair = proposeRandomNni(tree);
			performNni(tree, nniPair);
			return nniPair;
		}

		// Preparing species edge wise node counts
		double[] sp_counts = new double[recon.host.nnodes];
		for (Node node : recon.guest.nodes)
			sp_counts[recon.nodeMap[node.id]]++;

		// Taking Square of the number of edges per host tree edge
		for (int i = 0; i < sp_counts.length; i++)
			sp_counts[i] *= sp_counts[i];

		// Picking most informative species node
		double[] sp_count_cumulative = Mathematics.prepareCumulative(sp_counts);
		int sp_nodeId = 0;
		double randv = prng.nextDouble();
		while (randv > sp_count_cumulative[sp_nodeId]) {
			++sp_nodeId;
		}
		Node snode = recon.host.nodes.get(sp_nodeId);

		// Picking guest nodes
		ArrayList<Node> nodeList = new ArrayList<Node>();
		for (Node node : tree.nodes) {
			if (recon.nodeMap[node.id] == snode.id)
				nodeList.add(node);
		}

		// Picking two nodes for branch swaping
		int rv;
		boolean iterate;
		Node node1, node2;
		Node node_a = null, node_b = null;

		int loop_counter = 0;
		do {

			if (++loop_counter > 20) {
				Pair<Node, Node> nniPair;
				nniPair = proposeRandomNni(tree);
				performNni(tree, nniPair);
				return nniPair;
			}

			iterate = false;

			rv = irand(nodeList.size(), prng);
			node_a = nodeList.get(rv);
			if (node_a.isRoot() || node_a.parent.isRoot()) {
				iterate = true;
				continue;
			}

			rv = irand(nodeList.size(), prng);
			node_b = nodeList.get(rv);
			if (node_b.isRoot() || node_b.parent.isRoot() || node_b.parent.id == node_a.parent.id) {
				iterate = true;
				continue;
			}

			if (!node_a.isLeaf()) { // Test if node_b is a descendent of node_a

				boolean under_a = false;
				for (Node ptr = node_b.parent; ptr != null; ptr = ptr.parent) {
					if (ptr.id == node_a.id) {
						under_a = true;
						break;
					}
				}
				if (under_a) {
					iterate = true;
					continue;
				}
			}

			if (!node_b.isLeaf()) { // Test if node_a is a descendent of node_b
				boolean under_b = false;
				for (Node ptr = node_a.parent; ptr != null; ptr = ptr.parent) {
					if (ptr.id == node_b.id) {
						under_b = true;
						break;
					}
				}
				if (under_b) {
					iterate = true;
					continue;
				}
			}

		} while (iterate);

		// Both are different cherries
		node1 = node_a.parent;
		node2 = node_b.parent;

		// find child indexes
		int a = (node1.children.get(0) == node_a) ? 0 : 1;
		assert (node1.children.get(a) == node_a);

		int b = (node2.children.get(0) == node_b) ? 0 : 1;
		assert (node2.children.get(b) == node_b);

		// swap parent pointers
		node_a.parent = node2;
		node_b.parent = node1;

		// swap child pointers
		node2.children.set(b, node_a);
		node1.children.set(a, node_b);

		Pair<Node, Node> nniPair = new Pair<Node, Node>(node_a, node_b);
		return nniPair;
	}

	private static int irand(int max, PRNG prng) {
		int i = (int) (prng.nextDouble() * max);
		return (i == max) ? max - 1 : i;
	}

	// =============================================================================
	// Subtree pruning and regrafting (SPR)
	// Derived from Rasmussen et. al. 2011
	/*
	 * a = subtree e = newpos
	 * 
	 * BEFORE .... f d / \ c e / \ ... a b ... ...
	 * 
	 * AFTER
	 * 
	 * f d / \ b c ... / \ a e ... ...
	 * 
	 * Requirements: 1. a (subtree) is not root or children of root 2. e (newpos) is
	 * not root, a, descendant of a, c (parent of a), or b (sibling of a) 3. tree is
	 * binary
	 * 
	 */
	public static void performSpr(Tree tree, Node subtree, Node newpos) {
		Node a = subtree;
		Node e = newpos;

		Node c = a.parent;
		Node f = c.parent;
		int bi = (c.children.get(0) == a) ? 1 : 0;
		Node b = c.children.get(bi);
		int ci = (f.children.get(0) == c) ? 0 : 1;
		Node d = e.parent;
		int ei = (d.children.get(0) == e) ? 0 : 1;

		d.children.set(ei, c);
		c.children.set(bi, e);
		f.children.set(ci, b);

		b.parent = f;
		c.parent = d;
		e.parent = c;

		// b->dist += c->dist;
		// e->dist /= 2.0;
		// c->dist = e->dist;

	}

	/*
	 * What if e == f (also equivalent to NNI) this is OK
	 * 
	 * BEFORE
	 * 
	 * d / \ e ... / \ c ... / \ a b ... ...
	 * 
	 * AFTER d / \ c / \ a e ... / \ b ... ...
	 * 
	 * What if d == f (also equivalent to NNI) this is OK
	 * 
	 * BEFORE
	 * 
	 * f / \ c e / \ ... a b ... ...
	 * 
	 * AFTER
	 * 
	 * f / \ b c ... / \ a e ... ...
	 */

	/*
	 * Requirements: 1. a (subtree) is not root or children of root 2. e (newpos) is
	 * not root, a, descendant of a, c (parent of a), or b (sibling of a) 3. tree is
	 * binary
	 */
	public static Pair<Node, Node> proposeRandomSpr(Tree tree) {

		Node subtree, newpos;

		assert (tree.nnodes >= 5);

		// find subtree (a) to cut off (any node that is not root or child of root)
		int choice;
		do {
			choice = Mathematics.irand(tree.nnodes);
		} while (tree.nodes.get(choice).parent == null || tree.nodes.get(choice).parent.parent == null);
		Node a = tree.nodes.get(choice);
		subtree = a;

		// find sibling (b) of a
		Node c = a.parent;
		int bi = (c.children.get(0) == a) ? 1 : 0;
		Node b = c.children.get(bi);

		// choose newpos (e)
		Node e = null;
		do {
			choice = Mathematics.irand(tree.nnodes);
			e = tree.nodes.get(choice);

			// test if e is a valid choice
			if (e.parent == null || e == a || e == c || e == b)
				continue;

			// also test if e is a descendent of a
			boolean under_a = false;
			for (Node ptr = e.parent; ptr != null; ptr = ptr.parent) {
				if (ptr == a) {
					under_a = true;
					break;
				}
			}

			if (under_a)
				continue;

			break;
		} while (true);
		newpos = e;

		Pair<Node, Node> sprPair = new Pair<Node, Node>(subtree, newpos);
		return sprPair;

	}

	public static boolean validSpr(Tree tree, Node subtree, Node newpos) {
		Node a = subtree;
		Node e = newpos;

		// find sibling (b) of a
		Node c = a.parent;
		int bi = (c.children.get(0) == a) ? 1 : 0;
		Node b = c.children.get(bi);

		// test if a is a valid choice
		if (a.parent == null || a.parent.parent == null)
			return false;

		// test if e is a valid choice
		if (e.parent == null || e == a || e == c || e == b) {
			// printf("ERROR: a=%d, e=%d, b=%d, c=%d\n", a->name, e->name,
			// b->name, c->name);
			return false;
		}

		// also test if e is a descendent of a
		for (Node ptr = e.parent; ptr != null; ptr = ptr.parent)
			if (ptr == a)
				return false;

		return true;
	}

	public static Pair<Node, Node> proposeRandomRerooting(Tree tree) {

		// find edges for rerooting
		int choice;
		do {
			choice = Mathematics.irand(tree.nnodes);
		} while (tree.nodes.get(choice).isRoot());

		Node node1 = tree.nodes.get(choice);
		Node node2 = tree.nodes.get(choice).parent;

		Pair<Node, Node> edge = new Pair<Node, Node>(node1, node2);
		return edge;
	}

	public static void performRerooting(Tree tree, Pair<Node, Node> edge) {

		Node newroot = null;
		Node node1 = edge.first;
		Node node2 = edge.second;

		// determine new root
		if (node1.parent == node2)
			newroot = node1;
		else if (node2.parent == node1)
			newroot = node2;
		else if (node1.parent == tree.root || node2.parent == tree.root)
			// do nothing
			return;

		tree.reroot(newroot, true);
	}

	// if(this.level == Parameters.GENE_LEVEL)
	// System.out.println(this.lastOperationType);

	/*
	 * //public Tree tree; //this.tree = this.treeCache; private boolean
	 * assertRevert() { int [] array1 = this.tree.parray; int [] array2 =
	 * this.treeCache.parray; for(int i=0 ; i<array1.length ; i++) if(array1[i] !=
	 * array2[i]) return false; return true; }
	 * 
	 * // Undo topology change by performing NNI & SPR again
	 * if(this.lastOperationType.equalsIgnoreCase("NNI")) revertNNI();
	 * 
	 * if(this.lastOperationType.equalsIgnoreCase("SPR")) revertSPR();
	 * 
	 * if(this.lastOperationType.equalsIgnoreCase("REROOTING")) revertREROOTING();
	 * 
	 * 
	 * private void revertNNI() { performNni(this.tree, nniPair); }
	 * 
	 * private void revertSPR() { performSpr(tree, this.sprPair.first,
	 * this.sprNodec); }
	 * 
	 * private void revertREROOTING() { performRerooting(this.tree,this.oldRoot); }
	 * 
	 * public void setTree(Tree tree) { this.tree = tree; }
	 */

}
