package seqEvolution;

import common.LogDouble;

import phylogeny.Tree;

public class RKSeqEvoModel {

	/*
	 * Procedure to native call interface 1) cd
	 * /Users/auwnm/Documents/Jworkspace/DomainDLRS/bin 2) javah -classpath
	 * /Users/auwnm/Documents/Jworkspace/DomainDLRS/bin -d
	 * /Users/auwnm/Documents/cworkspace/RasmussenSeqEvo -jni seqEvolution.RKSeqEvo
	 * 3)
	 * 
	 */

	public static float ratio = (float) 1.0;
	public static float[] bgfreq = { (float) 0.25, (float) 0.25, (float) 0.25, (float) 0.25 };

	public static LogDouble getSeqLikelihood(Tree tree, MSAData msa) {
		double loglk;

		// if(msa.leafOrderedSeqs == null)
		msa.setLeafOrderedSequences(tree);

		// if(tree.parray == null)
		tree.parray = tree.tree2ptree();

		loglk = RKSeqEvo.calcSeqLikelihood(tree.nnodes, tree.parray, tree.bl, msa.getNoOfSequences(),
				msa.leafOrderedSeqs, bgfreq, ratio);
		return new LogDouble(loglk, 1);
	}

}
