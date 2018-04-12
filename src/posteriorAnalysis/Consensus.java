package posteriorAnalysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import jebl.evolution.trees.*;
import jebl.evolution.io.*;

public class Consensus {

	/*
	 * JEBL by Andrew Rambaut et al.
	 * 
	 * Builds greedy consensus tree given a set of unrooted trees. Each edge in a
	 * tree "supports" a split, i.e. the partition of the taxa into two clades
	 * (which you get by deleting the edge). A Majority consensus tree is built by
	 * finding clades appearing in at least 50% of the trees or more. A Greedy
	 * consensus tree is a refinement of a Majority tree where the splits are sorted
	 * by amount of support, and are applied in order only if they are compatible
	 * (consistent) with the tree at that stage. A user supplied threshold gives a
	 * lower bound on amount of support. If set t0 50% the Greedy method reduces to
	 * the majority consensus tree. At 100% it reduces to the Strict consensus tree.
	 * 
	 * The implementation is relatively simple but tricky in parts. Each tree is
	 * scanned, and support for each split/clade is collected in one table. The
	 * clade is represented by a bitset, which always contains the (arbitrary) first
	 * node. The scan is made by going over all nodes ordered in such a way that the
	 * subtree of exactly one edge of the node has been completely scanned, so the
	 * node "knows" the set of tips of that subtree without needing to re-scan the
	 * tree. After all trees are scanned an initial consensus tree is constructed
	 * with one root and all tips as children. The split set is scanned in order of
	 * decreasing support, and each supported clade refines the tree by creating a
	 * new descendant for the node containing the clade and re-attaching the clade
	 * to that new node. This is done only if the split is compatible with the tree,
	 * i.e. only if the split is completely contained in a proper subset of
	 * descendants of one node. This process continues until only clades with
	 * support lower that the threshold are left. The length of the consensus tree
	 * branches is computed from the average over all trees containing the clade.
	 * The lengths of tip branches are computed by averaging over all trees. While
	 * the consensus tree is logically unrooted, we generate a rooted tree because
	 * we can store attributes such as support only for nodes. Author: Joseph Heled
	 * Note: code to use JEBL was derived from phyutility by Dunn et al. Refer to
	 * Paths: http://sourceforge.net/projects/jebl/
	 * http://jebl.sourceforge.net/doc/api/jebl/evolution/trees/
	 * GreedyUnrootedConsensusTreeBuilder.html https://code.google.com/p/phyutility/
	 * https://github.com/blackrim/phyutility
	 */

	public static String performConsensus(String filename, Double threshold) throws IOException, ImportException {

		String treeStr = null;

		NewickImporter ni = new NewickImporter(new FileReader(filename), true);
		ArrayList<Tree> trees = (ArrayList<Tree>) ni.importTrees();
		RootedTree[] tr = new RootedTree[trees.size()];
		for (int i = 0; i < tr.length; i++) {
			tr[i] = (RootedTree) trees.get(i);
		}
		GreedyRootedConsensusTreeBuilder tb = new GreedyRootedConsensusTreeBuilder(tr, threshold, "greedy", true);
		treeStr = Utils.toNewick(Utils.rootTheTree(tb.build()));
		return treeStr;
	}

	public static String getConsensusTree(ArrayList<String> nwkStr) {

		String treeStr = null;
		try {

			// Setting threshold for consensus (1.0 = strict, 0.5 = majrule, 0 = allcompat)
			double threshold = 0.5;

			// Create temporary file to list the trees
			String tempFileName = "trees_for_consensus";
			File tempFile = File.createTempFile(tempFileName, ".tmp");
			BufferedWriter bwOut = new BufferedWriter(new FileWriter(tempFile));
			for (String tree : nwkStr) {
				bwOut.write(tree + "\n");
			}
			bwOut.flush();

			// Creating consensus tree
			treeStr = performConsensus(tempFile.getAbsolutePath(), Double.valueOf(threshold));

			// Deleting temporary files
			tempFile.delete();

			// Closing BufferedWriter
			bwOut.close();

		} catch (FileNotFoundException e) {
		} catch (IOException e) {
		} catch (ImportException e) {
		}

		return treeStr;
	}

}

// String outfile
// FileWriter wr = new FileWriter(outfile);
// NexusExporter ne = new NexusExporter(wr);
// ne.exportTree(tb.build());
// NewickExporter ne = new NewickExporter(wr);
// ne.exportTree(tb.build());
// wr.close();
/*
 * ArrayList<String> toyTree = new ArrayList<String>(); toyTree.add(
 * "(((one, two),(three,four)),(five,six));" ); toyTree.add(
 * "(((one, two),(three,five)),(four,six));" ); toyTree.add(
 * "(((one, two),(three,five)),(four,six));" ); toyTree.add(
 * "(((one, two),(three,five)),(four,six));" ); toyTree.add(
 * "(((one, two),(three,four)),(five,six));" ); // ( ( (two,one),(five,three)
 * ),(four,six))
 */
/*
 * ArrayList<String> toyTree = new ArrayList<String>(); toyTree.add(
 * "(((one, two),(three,four)),(five,six));" ); toyTree.add(
 * "(((one, two),(three,five)),(four,six));" ); toyTree.add(
 * "(((one, two),(three,five)),(four,six));" ); toyTree.add(
 * "(((one, two),(three,five)),(four,six));" ); toyTree.add(
 * "(((one, two),(three,four)),(five,six));" ); toyTree.add(
 * "(((one, two),(three,four)),(five,six));" );
 * System.out.println(Consensus.getConsensusTree(toyTree));
 */
