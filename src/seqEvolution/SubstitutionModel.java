package seqEvolution;

import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import phylogeny.Node;
import phylogeny.Tree;
import common.LogDouble;
import common.Pair;

/**
 * Implements the standard (probabilistic) Markov model for the substitution
 * process of sequence evolution for any type of aligned sequence data, see
 * e.g., Felsenstein 1981.
 * <p/>
 * Sequence positions with identical state patterns over leaves are identified
 * and counted. Rate variation across sites over discrete classes, e.g., Yang
 * 1993, can be modelled.
 * <p/>
 * This model does not yet support partitions of data into user-defined
 * "independent" loci (domains, codon positions, etc.). NOTE: This class is
 * derived from the C++ class <code>CacheSubstitutionModel</code> and not
 * <code>FastCacheSubstitutionModel</code> since the latter was stated
 * unsuitable for tree topology changes in the C++ CMake default settings. /Joel
 * 
 * @author Bengt Sennblad.
 * @author Lars Arvestad.
 * @author Joel Sjöstrand.
 */
public class SubstitutionModel {

	private MSAData D; // Sequence data (MSA).
	private SubstitutionMatrix Q; // Transition rate matrix Q and, with it, P
	private Tree T; // Rooted Binary Tree
	private LinkedHashMap<Integer, Double> L; // Associated lengths indexed w.r.t. node id;
	private boolean useRootEdge; // Decides if root arc should be included in computations.
	private int alphabetSize;

	/**
	 * For each vertex n of V(T), holds the likelihoods for the planted subtree T^n.
	 * Each such element is an array with rows corresponding to unique patterns and
	 * columns corresponding to site rate categories. Each element in such an array
	 * holds a vector with likelihoods corresponding to the states of the sequence
	 * type alphabet.
	 */

	public ArrayList<LikelihoodVectors> mylikelihoods;
	private LogDouble modelLikelihood; // Model likelihood.
	private DenseMatrix64F tmp; // Temporary matrix used during computations.

	/**
	 * Constructor.
	 * 
	 * @param D
	 *            sequence data (MSA).
	 * @param Q
	 *            data transition matrix Q (and P).
	 * @param T
	 *            tree.
	 * @param names
	 *            leaf names of T.
	 * @param branchLengths
	 *            branch lengths of T.
	 * @param useRootArc
	 *            if true, utilises the root arc ("stem") branch length when
	 *            computing model likelihood; if false, discards the root arc.
	 */
	public SubstitutionModel(MSAData D, SubstitutionMatrix Q, Pair<Tree, LinkedHashMap<Integer, Double>> T,
			boolean useRootEdge) {

		this.D = D;
		this.Q = Q;
		this.T = T.first;
		this.L = T.second;
		int noOfVertices = this.T.nnodes;
		this.useRootEdge = useRootEdge;
		this.alphabetSize = Q.getAlphabetSize();
		this.modelLikelihood = new LogDouble(0.0);
		this.tmp = new DenseMatrix64F(alphabetSize, 1);
		this.mylikelihoods = new ArrayList<LikelihoodVectors>(noOfVertices);
		for (int i = 0; i < noOfVertices; ++i) {
			LikelihoodVectors pv = new LikelihoodVectors(this.D.getNoOfPositions(), alphabetSize);
			mylikelihoods.add(pv);
		}

	}

	public SubstitutionModel(MSAData D, SubstitutionMatrix Q, Tree tree, boolean useRootEdge) {

		this.D = D;
		this.Q = Q;
		this.T = tree;

		LinkedHashMap<Integer, Double> lengths = new LinkedHashMap<Integer, Double>();
		for (int i = 0; i < tree.nnodes; i++)
			lengths.put(tree.nodes.get(i).id, tree.bl[tree.nodes.get(i).id]);

		this.L = lengths;
		int noOfVertices = this.T.nnodes;
		this.useRootEdge = useRootEdge;
		this.alphabetSize = Q.getAlphabetSize();
		this.modelLikelihood = new LogDouble(0.0);
		this.tmp = new DenseMatrix64F(alphabetSize, 1);
		this.mylikelihoods = new ArrayList<LikelihoodVectors>(noOfVertices);
		for (int i = 0; i < noOfVertices; ++i) {
			LikelihoodVectors pv = new LikelihoodVectors(this.D.getNoOfPositions(), alphabetSize);
			mylikelihoods.add(pv);
		}

	}

	public void updateTreeParameter(Tree newTree) {
		this.T = newTree;
	}

	private void computeModelLikelihood() {

		// Get root likelihood and patterns.
		Node node = this.T.root;

		// The likelihood data structures must be up-to-date.
		this.updateLikelihood(node, true);
		LikelihoodVectors pl = this.mylikelihoods.get(node.id);

		// Reset model likelihood.
		this.modelLikelihood = new LogDouble(1.0);

		// For every column in msa
		for (int j = 0; j < this.D.getNoOfPositions(); j++) {

			LogDouble patternL = new LogDouble(0.0);
			DenseMatrix64F curr = pl.likelihoods[j];

			// Multiply with stationary frequencies (that's our assumption for evolution
			// start).
			this.Q.multiplyWithPi(curr, this.tmp);
			patternL.add(CommonOps.elementSum(this.tmp));

			// Multiply with overall likelihood.
			this.modelLikelihood.mult(patternL);
		}
	}

	/**
	 * DP method which updates the likelihood column vectors for a subtree.
	 * 
	 * @param n
	 *            vertex root of subtree.
	 * @param doRecurse
	 *            true to process entire subtree rooted at n; false to only do n.
	 */
	private void updateLikelihood(Node node, boolean doRecurse) {
		if (node.isLeaf()) {
			this.updateLeafLikelihood(node);
		} else {

			// Process kids first.
			if (doRecurse) {
				this.updateLikelihood(node.children.get(0), true);
				this.updateLikelihood(node.children.get(1), true);
			}

			// Get data and likelihood storage.
			LikelihoodVectors pl = this.mylikelihoods.get(node.id);

			// Get child likelihoods.
			LikelihoodVectors pl_l = this.mylikelihoods.get(node.children.get(0).id);
			LikelihoodVectors pl_r = this.mylikelihoods.get(node.children.get(1).id);

			// Just a special case: we discard evolution over the stem arc if desired (when
			// doUseP = false).
			boolean doUseP = (this.useRootEdge || !node.isRoot());

			// Compute Pr[Dk | T, l, r(j)]

			if (doUseP) {
				// Set up site rate-specific P matrix.
				double w = L.get(node.id);
				this.Q.updateTransitionMatrix(w);
			}

			// Lastly, loop over each unique pattern in patterns.
			for (int i = 0; i < this.D.getNoOfPositions(); i++) {
				DenseMatrix64F left = pl_l.likelihoods[i];
				DenseMatrix64F right = pl_r.likelihoods[i];

				// Element-wise multiplication, tmp = left .* right.
				CommonOps.elementMult(left, right, this.tmp);
				DenseMatrix64F curr = pl.likelihoods[i];
				if (doUseP) {
					this.Q.multiplyWithP(this.tmp, curr);
				} else {
					pl.likelihoods[i].set(this.tmp);
				}
			}
		}

	}

	// Dynamic Programming Method which updates the likelihoods column vector for a
	// leaf vertex.
	private void updateLeafLikelihood(Node leaf) {

		// Set up data and likelihood storage.
		LikelihoodVectors pl = this.mylikelihoods.get(leaf.id);

		// Get sequence index for this vertex.
		int seqIdx = this.D.getSequenceIndex(leaf.name);

		// Set up site rate-specific P matrix.
		double w = L.get(leaf.id);
		this.Q.updateTransitionMatrix(w);

		// Loop over each unique pattern in patterns.
		for (int j = 0; j < this.D.getNoOfPositions(); j++) {
			// Compute likelihood.
			DenseMatrix64F curr = pl.likelihoods[j];
			int state = this.D.getIntState(seqIdx, j);
			this.Q.getLeafLikelihood(state, curr);
		}
	}

	public LogDouble getLogLikelihood() {
		this.computeModelLikelihood();
		return this.modelLikelihood;
	}

	public void printLogLikelihood() {
		this.computeModelLikelihood();
		System.out.println("P(D|T) = " + this.modelLikelihood.getLogValue());
	}

}
