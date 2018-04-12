package common;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;

import phylogeny.Node;
//import phylogeny.Reconciliation;
import phylogeny.Reconciliation;
import phylogeny.Tree;
import phylogeny.visTree;
import proposers.Proposer;
import proposers.ProposerSelector;
import proposers.ProposerStatistics;
import Main.Input;
import Main.Parameters;
import Main.mainSystem;

public class Log {

	private File logFile;
	private BufferedWriter bwLog;

	public Log(String fileName) throws IOException {
		if (fileName == null)
			System.out.println("Error: Log file name is null.");
		logFile = new File(fileName);
		bwLog = new BufferedWriter(new FileWriter(logFile));
	}

	public Log(String fileName, boolean append) throws IOException {
		logFile = new File(fileName);
		if (append && !logFile.exists()) {
			System.out.println("Error: Cannot read log file " + fileName);
		}
		bwLog = new BufferedWriter(new FileWriter(logFile));
	}

	public void commit() throws IOException {
		this.bwLog.flush();

	}

	public void write(String str) throws IOException {
		this.bwLog.write(str);

	}

	public void newLine() throws IOException {
		this.bwLog.newLine();
	}

	public void writePreRunInfo(Input input) throws IOException {

		this.write("# =========================================================================\n");
		this.write("# DomainDLRS Pre-Run Information                                           \n");
		this.write("# =========================================================================\n");
		Calendar cal = Calendar.getInstance();
		SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd hh:mm:ss");
		this.newLine();
		this.write("Current time : " + df.format(cal.getTime()) + '\n');
		this.newLine();
		this.write("Input arguments: \n");
		this.write("Number of domains  : " + input.numberOfDomains + "\n");
		this.write("Substitution Model : " + input.substitutionModel + "\n");
		this.write("Working Directory  : " + input.path + "\n");
		this.write("Species Tree File  : " + input.speciesTreeFile + "\n");
		this.write("Gene Map File      : " + input.geneMapFile + "\n");
		this.write("Gene Tree File     : " + input.geneTree + "\n");
		this.write("Output Prefix      : " + input.outPrefix + "\n");

		for (int i = 0; i < input.numberOfDomains; i++) {
			this.write("Domain msa file    : " + input.domAlignmentFiles[i] + "\n");
			this.write("Domain map file    : " + input.domMappingFiles[i] + "\n");
		}
	}

	// This method will write proposers statisitcs & acceptence count
	public void writePostRunInfo(ProposerSelector ps, ProposerStatistics runStats) throws IOException {

		this.newLine();
		this.write("# =========================================================================\n");
		this.write("# DomainDLRS Post-Run Information                                          \n");
		this.write("# =========================================================================\n");
		this.newLine();

		this.write("\nProposers statistics : \n");

		for (Proposer prop : ps.geneLevelProposers) {
			this.newLine();
			prop.getProposerStatistics().writeProposerStats(this);
			this.write("---------------------------------------\n");
		}
		for (int r = 0; r < ps.numberOfDomains; r++) {
			for (Proposer prop : ps.domLevelProposers.get(r)) {
				this.newLine();
				prop.getProposerStatistics().writeProposerStats(this);
				this.write("---------------------------------------\n");
			}
		}
		this.newLine();

		this.write("\nRun statistics : \n");
		runStats.writeProposerStats(this);

	}

	// This method will write proposers statisitcs & acceptence count
	public void writeBestState(StringBuilder bestSample) throws IOException {
		this.newLine();
		this.write("# =========================================================================\n");
		this.write("# Best encountered state:                                                  \n");
		this.write("# =========================================================================\n");
		this.newLine();
		this.write(bestSample.toString());
	}

	// This method will log the tree parameters i.e. species,gene and domain trees
	public void writeTreeParameters(Parameters pi, boolean drawTrees) throws IOException {

		writeTree(pi.speciesTree, true, null, null);
		writeTree(pi.geneTree, true, null, null);
		for (int i = 0; i < pi.numberOfDomains; i++)
			writeTree(pi.domainTree[i], true, "len", pi.domainTree[i].bl);

		if (drawTrees) {
			visTree.writeTree(pi.speciesTree, mainSystem.CWD + "speciesTree");
			visTree.writeTree(pi.geneTree, mainSystem.CWD + "geneTree");

			for (int i = 0; i < pi.numberOfDomains; i++)
				visTree.writeTree(pi.domainTree[i], mainSystem.CWD + "domTree" + (i + 1));
		}
	}

	// This method will write the rates parameters i.e. Birth death and edge rates
	public void writeRateParameters(Parameters pi) throws IOException {

		String rateStr;
		StringBuilder sb = new StringBuilder();

		rateStr = String.format("[%.5f , %.5f]", pi.geneBDRates[0], pi.geneBDRates[1]);
		sb.append("\nGene BD rates = " + rateStr);
		for (int i = 0; i < pi.numberOfDomains; i++) {
			rateStr = String.format("[%.5f , %.5f]", pi.domBDRates[i][0], pi.domBDRates[i][1]);
			sb.append("\nDomain" + (i + 1) + " BD rates = " + rateStr);
		}

		sb.append("\n");
		for (int i = 0; i < pi.numberOfDomains; i++) {
			rateStr = String.format("[%.5f , %.5f]", pi.domEdgeRates[i][0], pi.domEdgeRates[i][1]);
			sb.append("\nDomain" + (i + 1) + " edge rate parameters = " + rateStr);
		}

		sb.append("\n");
		bwLog.write(sb.toString());
	}

	// This method will write the reconciliation between species, gene
	public void writeReconciliations(Reconciliation geneRecon, Reconciliation[] domainRecon) throws IOException {
		writeReconciliation(geneRecon);
		for (int i = 0; i < domainRecon.length; i++)
			writeReconciliation(domainRecon[i]);
	}

	public void writeParameters(Parameters pi, Reconciliation geneRecon, Reconciliation[] domainRecon,
			boolean drawTrees) throws IOException {

		String rateStr;
		StringBuilder sb = new StringBuilder();

		rateStr = String.format("[%.5f , %.5f]", pi.geneBDRates[0], pi.geneBDRates[1]);
		sb.append("\nGene BD rates = " + rateStr);
		for (int i = 0; i < pi.numberOfDomains; i++) {
			rateStr = String.format("[%.5f , %.5f]", pi.domBDRates[i][0], pi.domBDRates[i][1]);
			sb.append("\nDomain" + (i + 1) + " BD rates = " + rateStr);
		}

		sb.append("\n");
		for (int i = 0; i < pi.numberOfDomains; i++) {
			rateStr = String.format("[%.5f , %.5f]", pi.domEdgeRates[i][0], pi.domEdgeRates[i][1]);
			sb.append("\nDomain" + (i + 1) + " edge rate parameters = " + rateStr);
		}

		sb.append("\n");
		bwLog.write(sb.toString());

		writeTree(pi.speciesTree, true, null, null);
		writeTree(pi.geneTree, true, null, null);
		for (int i = 0; i < pi.numberOfDomains; i++)
			writeTree(pi.domainTree[i], true, "len", pi.domainTree[i].bl);

		if (geneRecon != null)
			writeReconciliation(geneRecon);

		if (domainRecon != null) {
			for (int i = 0; i < pi.numberOfDomains; i++)
				writeReconciliation(domainRecon[i]);
		}
		if (drawTrees) {
			visTree.writeTree(pi.speciesTree, mainSystem.CWD + "speciesTree");
			visTree.writeTree(pi.geneTree, mainSystem.CWD + "geneTree");

			for (int i = 0; i < pi.numberOfDomains; i++)
				visTree.writeTree(pi.domainTree[i], mainSystem.CWD + "domTree" + (i + 1));
		}
	}

	// This method will write different measurement w.r.t vertices
	// e.g. branch Length, probabilities etc.
	public void writeTree(Tree tree, boolean postOrder, String measureName, double[] measure) throws IOException {

		StringBuilder sb = new StringBuilder();

		String child1 = new String();
		String child2 = new String();

		sb.append("\n" + tree.name + "\n");
		sb.append("------------------------\n");

		// Please check if i did it correctly ....

		ArrayList<Node> nodeList;
		nodeList = tree.nodes;

		if (postOrder) {
			nodeList = new ArrayList<Node>();
			tree.getTreePostOrder(nodeList, tree.root);
		}

		for (Node node : nodeList) {
			if (node.isLeaf())
				sb.append("Leaf " + String.format("%3d", node.id) + " " + String.format("%-30s", node.name));
			else {
				sb.append("Node " + String.format("%3d", node.id) + " ");
				for (int j = 0; j < node.nchildren; j++) {
					if (j != node.nchildren - 1)
						child1 = String.format("%-3d", node.children.get(j).id);
					else
						child2 = String.format("%-3d", node.children.get(j).id);
				}
				String interName = "(" + child1 + "," + child2 + ")";
				sb.append(String.format("%-30s", interName));
			}

			if (tree.vt != null) {
				sb.append(" vt = " + String.format("%3.8f", tree.vt[node.id]));
			}
			if (measure != null) {
				sb.append(measureName + " = " + String.format("%3.8f", measure[node.id]));
			}
			if (node.isRoot())
				sb.append("\t root node");

			sb.append("\n");
		}
		this.write(sb.toString());

	}

	// Helping function for writing log of a single iteration probabilities and sum
	public void writeHillIteration(int nodeid, double position, double delta, LogDouble currentLikelihood,
			LogDouble previousLikelihood, boolean header) throws IOException {

		if (header) {
			String str = String.format("%15s", "nodeid");

			str += String.format("\t%-15s", "position");

			str += String.format("\t%-15s", "Change");

			str += String.format("\t%-20s", "Current Likelihood");

			str += String.format("\t%-20s", "Previous Likelihood");

			this.write(str);
			this.newLine();
		}

		String likelihoodsStr = String.format("%5d\t", nodeid);

		likelihoodsStr += String.format("%3.8f\t", position);

		likelihoodsStr += String.format("%3.8f\t", delta);

		likelihoodsStr += String.format("%3.8f\t", currentLikelihood.getLogValue());

		likelihoodsStr += String.format("%3.8f\t", previousLikelihood.getLogValue());

		this.write(likelihoodsStr);
		this.newLine();
	}

	// Helping function for writing log of a single iteration probabilities and sum
	public void writeIteration(int iteration, Parameters pi, LogDouble geneDLModelLikelihood,
			LogDouble[] domainDLModelLikelihood, LogDouble[] rateModelLikelihoods, LogDouble Product)
			throws IOException {

		if (iteration == 0) {
			String str = String.format("%5s", "i");

			str += String.format("\t%-15s", "P(T_G|theta)");

			for (int i = 0; i < pi.numberOfDomains; i++)
				str += String.format("\t%-15s", "P(T_D" + (i + 1) + "|theta)");

			for (int i = 0; i < pi.numberOfDomains; i++)
				str += String.format("\t%-15s", "P(l|t_D,T_D" + (i + 1) + ")");

			str += String.format("\t%-15s", "Product");

			this.write(str);
			this.newLine();
		}

		String likelihoodsStr = String.format("%5d\t", iteration);
		likelihoodsStr += String.format("%3.8f\t", geneDLModelLikelihood.getLogValue());
		for (int i = 0; i < pi.numberOfDomains; i++)
			likelihoodsStr += String.format("%3.8f\t", domainDLModelLikelihood[i].getLogValue());
		for (int i = 0; i < pi.numberOfDomains; i++)
			likelihoodsStr += String.format("%e\t", rateModelLikelihoods[i].getLogValue());

		likelihoodsStr += String.format("%e\t", Product.getLogValue());

		this.write(likelihoodsStr);
		this.newLine();
	}

	public boolean writeExtinctionProb(Tree htree, double[] extinctionPrb) throws IOException {
		String str;
		StringBuilder sb = new StringBuilder();
		ArrayList<Node> orderedNodes = new ArrayList<Node>();
		htree.getTreePostOrder(orderedNodes, htree.root);

		sb.append("Death Probabilities along " + htree.name + ":\n");
		for (int i = 0; i < htree.nnodes; i++) {
			if (orderedNodes.get(i).isLeaf())
				str = String.format("@ Leaf %-3s %-25s", orderedNodes.get(i).id, orderedNodes.get(i).name);
			else
				str = String.format("@ Node %-3s parent(%-3s,%-3s)\t", orderedNodes.get(i).id,
						orderedNodes.get(i).children.get(0).id, orderedNodes.get(i).children.get(1).id);
			sb.append(str + "\tP(Extinction) = " + Math.exp(extinctionPrb[orderedNodes.get(i).id]) + "\n");
		}

		this.newLine();
		this.write(sb.toString());
		return true;
	}

	public boolean writeReconciliation(Reconciliation R) throws IOException {

		StringBuilder sb = new StringBuilder();
		sb.append("MPR Reconciliation between " + R.guest.name + " and " + R.host.name + ":\n");
		int nnodes = R.guest.nnodes;
		for (int i = 0; i < nnodes; i++) {
			String guestID = String.format("%-2d", R.guest.nodes.get(i).id);
			String hostID = String.format("%-2d", R.nodeMap[R.guest.nodes.get(i).id]);
			int event = R.eventMap[R.guest.nodes.get(i).id];

			String showEvent = null;
			if (event == Reconciliation.EVENT_LEAF)
				showEvent = "Leaf";
			else if (event == Reconciliation.EVENT_SPEC)
				showEvent = "Speciation";
			else if (event == Reconciliation.EVENT_DUPL)
				showEvent = "Duplication";

			if (showEvent != null)
				sb.append(guestID + " --> " + hostID + "\t" + showEvent + "\n");
			else
				sb.append(guestID + " --> " + hostID + "\n");
		}

		this.newLine();
		this.write(sb.toString());
		return true;
	}

	public void close() throws IOException {
		bwLog.close();
	}

}
