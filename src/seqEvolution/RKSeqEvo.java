package seqEvolution;

public class RKSeqEvo {

	public static native boolean isThere();

	public static native int[] testParameter(int[] param);

	public static native double calcSeqLikelihood(int nnodes, int[] ptree, double[] dists, int nseqs, String[] seqs,
			float[] bgfreq, float ratio);

}
