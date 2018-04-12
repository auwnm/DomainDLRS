package seqEvolution;

import java.util.HashMap;
import org.ejml.alg.dense.decomposition.DecompositionFactory;
import org.ejml.alg.dense.decomposition.EigenDecomposition;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import common.LogDouble;

/**
 * Handles transition probabilities of a Markov process for molecular sequence
 * evolution. As such, it can be viewed as a representing a momentary transition
 * rate matrix Q, as well as a transition probability matrix P for some user-set
 * Markov time w.
 * <p/>
 * The momentary transition rate matrix Q can be decomposed into a symmetric
 * <i>exchangeability</i> matrix R and the vector of stationary frequencies Pi.
 * R is symmetric 'intrinsic rate' matrix of the model (a.k.a. 'exchangeability'
 * matrix), Represented as a row-major triangular matrix with size
 * (dim*(dim-1)/2,1). Symmetric values are implicitly defined due to time
 * reversibility assumption. The transition probability matrix P over a given
 * (Markov) time interval w is given by P=exp(Qw). Note that w often is measured
 * in the expected number of events per site occurring over the interval.
 * <p/>
 * Assumes time reversibility, that is, pi_i*mu_ij=pi_j*mu_ji for stationary
 * frequencies pi_x and transition probabilities mu_xy.
 * 
 * @author Bengt Sennblad.
 * @author Lars Arvestad.
 * @author Joel Sjöstrand.
 */
public class SubstitutionMatrix {

	public static final double MAX_MARKOV_TIME = 1000.0; // The maximum allowed time w on which transition rate matrix Q
															// can act.
	private String modelName; // Substitution model name.
	private SequenceType sequenceType; // The sequence type that is handled, e.g., DNA. TODO: Bens: Can we avoid it? */
	private int alphabetSize; // Number of states in model alphabet, dim.
	private DenseMatrix64F R; // The symmetric 'intrinsic rate' matrix of the model (a.k.a. 'exchangeability'
								// matrix)
	private DenseMatrix64F Pi; // The stationary frequencies of the model formulated as a diagonal matrix. Only
								// diagonal stored as matrix of size (dim,1).
	private DenseMatrix64F Q; // The transition rate matrix Q, normalised to have 1 expected event over branch
								// length 1.
	private DenseMatrix64F E; // Eigenvalues of the transition rate matrix. Only diagonal, stored as matrix of
								// size (dim,1).
	private DenseMatrix64F V; // Eigenvectors of the transition rate matrix. Size (dim,dim).
	private DenseMatrix64F iV; // Inverse of V. Size (dim,dim).
	private DenseMatrix64F P; // Transition probability matrix (updated frequently for different user times
								// w). Size (dim,dim).
	private DenseMatrix64F dP; // Derivative of transition probability matrix P(w) w.r.t. time variable w.
	private DenseMatrix64F tmp_matrix; // Temporary storage matrix. Size (dim,dim).
	private DenseMatrix64F tmp_diagonal; // Temporary storage vector. Size (dim,1).

	/**
	 * Small cache for ambiguity leaf likelihoods. Cleared every time P i updated.
	 */
	private HashMap<Integer, DenseMatrix64F> ambigCache;

	/**
	 * Constructor.
	 * 
	 * @param modelName
	 *            name of substitution model.
	 * @param sequenceType
	 *            sequence type.
	 * @param R_vec
	 *            'intrinsic' rate matrix, Row-major format, Time reversibility,
	 *            Length dim*(dim-1)/2, dim= length(alphabet)
	 * @param Pi_vec
	 *            stationary frequencies. Should have length dim, where dim is the
	 *            alphabet length.
	 * @param cacheSize
	 *            number of P matrices to store in cache, e.g., 1000.
	 */
	public SubstitutionMatrix(String modelName, SequenceType sequenceType, double[] R_vec, double[] Pi_vec) {
		this.modelName = modelName;
		this.sequenceType = sequenceType;
		this.alphabetSize = sequenceType.getAlphabetSize();
		this.R = new DenseMatrix64F(alphabetSize * (alphabetSize - 1) / 2, 1, true, R_vec);
		this.Pi = new DenseMatrix64F(alphabetSize, 1, true, Pi_vec);
		this.Q = new DenseMatrix64F(alphabetSize, alphabetSize);
		this.E = new DenseMatrix64F(alphabetSize, 1);
		this.V = new DenseMatrix64F(alphabetSize, alphabetSize);
		this.iV = new DenseMatrix64F(alphabetSize, alphabetSize);
		this.P = new DenseMatrix64F(alphabetSize, alphabetSize);
		this.tmp_matrix = new DenseMatrix64F(alphabetSize, alphabetSize);
		this.tmp_diagonal = new DenseMatrix64F(alphabetSize, 1);
		this.ambigCache = new HashMap<Integer, DenseMatrix64F>(14); // Not more than at most ~14 different ambiguity
																	// characters.
		this.update();
	}

	// Returns the name of the substitution model.
	public String getModel() {
		return this.modelName;
	}

	// Returns the sequence type that the model handles (DNA, AA, etc.).
	public SequenceType getSequenceType() {
		return this.sequenceType;
	}

	// Returns the alphabet size of the Markov process.
	public int getAlphabetSize() {
		return this.alphabetSize;
	}

	// Updates Q and the eigensystem based on R and Pi.
	private void update() {
		// Creates Q by means of R and Pi.
		// The diagonal values of Q = -the sum of other values of row, by definition.
		// R in this implementation holds upper triangle of symmetric matrix, excluding
		// diagonal.
		this.ambigCache.clear();
		this.Q.zero();
		int R_i = 0;
		double val;
		for (int i = 0; i < alphabetSize; i++) {
			for (int j = i + 1; j < alphabetSize; j++) {
				val = this.Pi.get(i) * this.R.get(R_i);
				this.Q.set(i, j, val);
				this.Q.set(i, i, this.Q.get(i, i) - val);
				// R is symmetric.
				val = this.Pi.get(j) * this.R.get(R_i++);
				this.Q.set(j, i, val);
				this.Q.set(j, j, this.Q.get(j, j) - val);
			}
		}

		// Perform scaling of Q so that a branch length of w=1 yields 1 expected event.
		double beta = 0;
		for (int i = 0; i < this.alphabetSize; ++i) {
			beta -= this.Pi.get(i) * this.Q.get(i, i);
		}
		beta = 1.0 / beta;
		CommonOps.scale(beta, this.Q);

		// Eigendecomposition of a matrix Q
		// NOTE: It is assumed solutions with imaginary parts will never be encountered.
		// To avoid checks, we assume non-symmetric Q. Symmetric models (JC69, etc.) are
		// rarely used in practice anyway.
		// EigenDecomposition<DenseMatrix64F> eigFact =
		// DecompositionFactory.eigSymm(this.alphabetSize, true);
		EigenDecomposition<DenseMatrix64F> eigFact = DecompositionFactory.eigGeneral(this.alphabetSize, true);
		if (!eigFact.decompose(this.Q)) {
			throw new RuntimeException("Unable to decompose eigensystem for substitution model.");
		}
		AdditionalEJMLOps.getEigensystemSolution(this.alphabetSize, eigFact, this.E, this.V);

		// Compute inverse of V.
		if (!CommonOps.invert(this.V, this.iV)) {
			throw new RuntimeException("Unable to invert matrix of eigenvectors in substitution model.");
		}
	}

	/**
	 * Sets up P=exp(Qt), the transition probability matrix for the Markov process
	 * over 'time' t (where 'time' is not necessarily temporal). Precondition: t <=
	 * 1000.
	 * 
	 * @param t
	 *            the "time" (or branch length) over which Q acts.
	 */
	public void updateTransitionMatrix(double t) {
		// C++ comment which may still apply... /joelgs
		// If w is too big, the precision of LAPACK seem to get warped!
		// The choice of max value of 1000 is arbitrary and well below the
		// actual max value. /bens
		// TODO: Could we precondition on a reasonable MarkovTime?
		if (t > MAX_MARKOV_TIME) {
			throw new IllegalArgumentException(
					"In substitution model, cannot compute transition probability matrix P for too large Markov time w="
							+ t + ".");
		}

		// Clear ambiguity cache.
		this.ambigCache.clear();
		this.tmp_diagonal.zero();
		this.tmp_matrix.zero();

		// P(t) = e^{Qt}
		AdditionalEJMLOps.elementExp(this.alphabetSize, this.E, t, this.tmp_diagonal);
		AdditionalEJMLOps.multDiagA(this.alphabetSize, this.tmp_diagonal, this.iV, this.tmp_matrix);
		this.P = new DenseMatrix64F(this.alphabetSize, this.alphabetSize);
		CommonOps.mult(this.V, this.tmp_matrix, this.P);

		// P(t) = e^{Qt}
		// Derivative: P'(t) = Q e^{Qt}
		this.dP = new DenseMatrix64F(this.alphabetSize, this.alphabetSize);
		CommonOps.mult(this.Q, this.P, this.dP);
	}

	public void displayProbTransitionMatrix() {

		for (int i = 0; i < this.P.numRows; i++) {
			for (int j = 0; j < this.P.numCols; j++)
				System.out.print(this.P.get(i, j) + "\t");
			System.out.println();
		}

	}

	/*
	 * Derivative of the loglikelihood (Adapted from FastPhylo tool) log L(t) = log
	 * \prod_{i,j} p_{i,j}^{N_{i,j}}(t) = \sum N_{i,j} log(p_{i,j}(t)) (log L(t))' =
	 * \sum_{i,j} N_{i,j} \cdot \frac{1}{p_{i,j}(t)} \cdot p_{i,j}'t() Which means
	 * that it becomes log L(t) = \sum N .* P'(t) ./ P(t) where P(t) = e^{Qt} P'(t)
	 * = Q e^{Qt}
	 */

	// This method will compute the value of derivative of log likelihood function
	// for a given replacement count matrix. i.e. log L(t) = \sum N .* P'(t) ./ P(t)
	public double getSeqDerivLogLikelihood(double t, double N[][]) {

		// Initializing variable for the function
		double likelihoodVal = 0;

		// Update the probability transition matrix P(t) & its derivative P'(t)
		updateTransitionMatrix(t);

		// Loop over each position in this alignment
		for (int i = 0; i < this.alphabetSize; i++)
			for (int j = 0; j < this.alphabetSize; j++)
				likelihoodVal += (N[i][j] * dP.get(i, j) / P.get(i, j));

		return likelihoodVal;
	}

	public double getSeqLogLikelihood(double t, double N[][]) {

		// Initializing likelihood for these two sequences.
		LogDouble likelihoodVal = new LogDouble(1.0);

		// Update the probability transition matrix P(t) & its derivative P'(t)
		updateTransitionMatrix(t);

		// Loop over each position in this alignment
		for (int i = 0; i < this.alphabetSize; i++)
			for (int j = 0; j < this.alphabetSize; j++) {
				LogDouble temp = new LogDouble(P.get(i, j));
				likelihoodVal.mult(N[i][j] * temp.getLogValue());
			}
		return likelihoodVal.getLogValue();
	}

	/**
	 * Performs matrix multiplication Y=P*X for the current P.
	 * 
	 * @param X
	 *            operand matrix (typically vector) of size (dim,ncol).
	 * @param Y
	 *            resulting matrix Y=P*X. Should have size (dim,ncol).
	 */
	public void multiplyWithP(DenseMatrix64F X, DenseMatrix64F Y) {
		CommonOps.mult(this.P, X, Y);
	}

	/**
	 * Returns the likelihood for a certain leaf state for the current P. This
	 * corresponds to the state's column in P (and analogously for ambiguity
	 * characters.
	 * 
	 * @param state
	 *            the state's integer index.
	 * @param result
	 *            the column values. Should have size (dim,1).
	 */
	public void getLeafLikelihood(int state, DenseMatrix64F result) {
		if (state < this.alphabetSize) {
			for (int i = 0; i < this.alphabetSize; ++i) {
				result.set(i, this.P.get(i, state));
			}
		} else {
			// Ambiguity state.
			DenseMatrix64F res = this.ambigCache.get(state);
			if (res == null) {
				// Not computed before.
				this.multiplyWithP(this.sequenceType.getLeafLikelihood(state), result);
				this.ambigCache.put(state, new DenseMatrix64F(result));
			} else {
				// Computed before.
				result.set(res);
			}
		}
	}

	/**
	 * Element-wise multiplication Y=Pi*X.
	 * 
	 * @param X
	 *            operand matrix (typically vector) of size (dim,ncol).
	 * @param Y
	 *            resulting matrix Y=Pi*X. Should have size (dim,ncol).
	 */
	public void multiplyWithPi(DenseMatrix64F X, DenseMatrix64F Y) {
		AdditionalEJMLOps.multDiagA(this.alphabetSize, this.Pi, X, Y);
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Markov transition rate matrix of model ").append(modelName).append('\n');
		sb.append("Current symmetric intrinsic rate matrix, R (flattened):\n");
		sb.append(this.R.toString()).append('\n');
		sb.append("Current stationary distribution base frequencies, Pi:\n");
		sb.append(this.Pi.toString()).append('\n');
		sb.append("Current eigenvalue matrix of Q, E:\n");
		sb.append(this.E.toString()).append('\n');
		sb.append("Current right-hand side eigenvectors of Q, V:\n");
		sb.append(this.V.toString()).append('\n');
		sb.append("Current inverse of V, iV:\n");
		sb.append(this.iV.toString()).append('\n');
		return sb.toString();
	};

	/**
	 * For debugging and similar purposes, returns R in a String format.
	 * 
	 * @return the R matrix.
	 */
	public String getRString() {
		StringBuilder sb = new StringBuilder();
		int R_index = 0;
		sb.append("Alphabet_size: ").append(alphabetSize).append('\n');
		for (int i = 0; i < alphabetSize; i++) {
			for (int j = 0; j < alphabetSize; j++) {
				if (j < alphabetSize) {
					sb.append('\t');
				}
				if (j > i) {
					sb.append(this.R.get(R_index++));
				}
			}
			if (i < alphabetSize - 2) {
				sb.append('\n');
			}
		}
		return sb.toString();
	}

}
