package treeConstructor;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map.Entry;

import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.evoinference.distance.NeighborJoining;
import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;

import common.PRNG;

import phylogeny.NewickReader;
import phylogeny.Tree;

public class NJConstructor {

	private static final DecimalFormat SAFE_NAME_NUMBER = new DecimalFormat("000000");

	// Tree gTree = NJConstructor.ForesterNJ (distances,
	// msa.getSequencesIdentifiers());
	// NJ tree construction by using forestor API
	public static Tree ForesterNJ(double[][] distances, String[] identifier) throws Exception {

		if (distances == null)
			System.out.println("Error: Distance matrix is not initialized");

		int IntID = 0;
		HashMap<String, String> idMap = new HashMap<String, String>();
		BasicSymmetricalDistanceMatrix matrix = new BasicSymmetricalDistanceMatrix(distances.length);
		;

		for (int i = 0; i < matrix.getSize(); i++) {
			String orgID = identifier[i];
			String newID = "NJID" + SAFE_NAME_NUMBER.format(IntID++);
			idMap.put(newID, orgID);
			matrix.setIdentifier(i, newID);
		}
		for (int col = 0; col < matrix.getSize(); col++) {
			for (int row = 0; row < matrix.getSize(); row++) {
				matrix.setValue(col, row, distances[col][row]);
			}
		}

		final NeighborJoining nj = NeighborJoining.createInstance();
		final Phylogeny p = nj.execute(matrix);

		boolean simpleNewick = true;
		boolean writeDistanceToParent = true;
		PhylogenyWriter pw = new PhylogenyWriter();
		String newickString = pw.toNewHampshire(p, simpleNewick, writeDistanceToParent).toString();

		for (Entry<String, String> name : idMap.entrySet()) {
			newickString = newickString.replaceAll(name.getKey(), name.getValue());
		}

		// Constructing tree object
		Tree tree = NewickReader.readNewickTreeStr(newickString.toString());
		return tree;
	}

	public static Tree RasmussenNJ(double distmat[][], String leavesName[]) {
		int ngenes = distmat.length;
		int nnodes = 2 * ngenes - 1;

		int[] ptree = new int[nnodes];
		double[] branches = new double[nnodes];
		double[][] dists = new double[nnodes][nnodes];
		double[] restdists = new double[nnodes];
		int[] leaves = new int[ngenes];
		int nleaves = ngenes;
		int newnode = ngenes;

		// initialize distances
		for (int i = 0; i < ngenes; i++) {
			double r = 0.0;
			for (int j = 0; j < ngenes; j++) {
				dists[i][j] = distmat[i][j];
				r += distmat[i][j];
			}
			restdists[i] = r / (ngenes - 2);
		}

		// initialize leaves
		for (int i = 0; i < ngenes; i++)
			leaves[i] = i;

		// join loop
		while (nleaves > 2) {
			// search for closest genes
			double low = Double.POSITIVE_INFINITY;
			int lowi = -1, lowj = -1;

			for (int i = 0; i < nleaves; i++) {
				for (int j = i + 1; j < nleaves; j++) {
					int gene1 = leaves[i];
					int gene2 = leaves[j];
					double dist = dists[gene1][gene2] - restdists[gene1] - restdists[gene2];
					if (dist < low) {
						low = dist;
						lowi = i;
						lowj = j;
					}
				}
			}

			// join gene1 and gene2
			int lowgene1 = leaves[lowi];
			int lowgene2 = leaves[lowj];
			int parent = newnode++;
			ptree[lowgene1] = parent;
			ptree[lowgene2] = parent;

			// set distances
			branches[lowgene1] = (dists[lowgene1][lowgene2] + restdists[lowgene1] - restdists[lowgene2]) / 2.0;
			branches[lowgene2] = dists[lowgene1][lowgene2] - branches[lowgene1];

			// gene1 and gene2 are no longer leaves, remove them from leaf set
			leaves[lowi] = parent;
			leaves[lowj] = leaves[nleaves - 1];
			nleaves--;

			double r = 0;
			for (int i = 0; i < nleaves; i++) {
				int gene = leaves[i];
				if (gene != parent) {
					double v = (dists[lowgene1][gene] + dists[lowgene2][gene] - dists[lowgene1][lowgene2]) / 2.0;
					dists[parent][gene] = v;
					dists[gene][parent] = v;
					r += v;
				}
			}

			if (nleaves > 2)
				restdists[parent] = r / (nleaves - 2);
		}

		// join the last two genes, split the remaining dist evenly
		int gene1 = leaves[0];
		int gene2 = leaves[1];
		int parent = newnode++;

		ptree[gene1] = parent;
		ptree[gene2] = parent;
		ptree[parent] = -1;
		branches[gene1] = dists[gene1][gene2] / 2.0;
		branches[gene2] = dists[gene1][gene2] / 2.0;
		branches[parent] = 0.0;

		assert (parent == ngenes * 2 - 2);

		// Constructing tree object
		Tree tree = Tree.ptree2tree(nnodes, ptree);
		tree.setLeafNames(leavesName, true);
		tree.setLengths(branches);

		/*
		 * LinkedHashMap<Integer,Double> lengths = new LinkedHashMap<Integer,Double>();
		 * for (int i=0; i<nnodes; i++) lengths.put(i,branches[i]); Pair<Tree ,
		 * LinkedHashMap<Integer,Double> > T = new Pair<Tree ,
		 * LinkedHashMap<Integer,Double> >(tree,lengths); Pair<Tree , double [] > T =
		 * new Pair<Tree , double [] >(tree,branches);
		 */

		return tree;
	}

	public static Tree RandomTree(String leavesName[], PRNG prng, double lengthScale) {
		int ngenes = leavesName.length;
		int nnodes = 2 * ngenes - 1;

		int[] ptree = new int[nnodes];
		double[] branches = new double[nnodes];
		int[] leaves = new int[ngenes];
		int nleaves = ngenes;
		int newnode = ngenes;

		// initialize leaves
		for (int i = 0; i < ngenes; i++)
			leaves[i] = i;

		// join loop
		while (nleaves > 2) {
			// search for pair of random genes
			int lowi = -1, lowj = -1;

			do {
				lowi = prng.nextInt(nleaves);
				lowj = prng.nextInt(nleaves);
			} while (lowi == lowj);

			// join gene1 and gene2
			int lowgene1 = leaves[lowi];
			int lowgene2 = leaves[lowj];
			int parent = newnode++;
			ptree[lowgene1] = parent;
			ptree[lowgene2] = parent;

			// set distances
			branches[lowgene1] = prng.nextDouble() * lengthScale;
			branches[lowgene2] = prng.nextDouble() * lengthScale;

			// gene1 and gene2 are no longer leaves, remove them from leaf set
			leaves[lowi] = parent;
			leaves[lowj] = leaves[nleaves - 1];
			nleaves--;

		}

		// join the last two genes, split the remaining dist evenly
		int gene1 = leaves[0];
		int gene2 = leaves[1];
		int parent = newnode++;

		ptree[gene1] = parent;
		ptree[gene2] = parent;
		ptree[parent] = -1;
		branches[gene1] = prng.nextDouble() * lengthScale;
		branches[gene2] = prng.nextDouble() * lengthScale;
		branches[parent] = prng.nextDouble() * lengthScale;

		assert (parent == ngenes * 2 - 2);

		// Constructing tree object
		Tree tree = Tree.ptree2tree(nnodes, ptree);
		tree.setLeafNames(leavesName, true);
		tree.setLengths(branches);

		return tree;
	}

}
