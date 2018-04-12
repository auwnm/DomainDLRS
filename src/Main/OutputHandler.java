package Main;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import mcmc.StateLikelihood;

import phylogeny.NewickWriter;
import phylogeny.Node;

import common.LogDouble;

public class OutputHandler {

	private BufferedWriter[] bwOut;

	public OutputHandler(int numberOfDomains, String path, String outPrefix) throws IOException {

		String fileName;
		File outputFile;
		bwOut = new BufferedWriter[numberOfDomains + 1];

		int r = 0;
		while (r < numberOfDomains) {
			fileName = path + outPrefix + "_dom" + r + ".mcmc";
			outputFile = new File(fileName);
			bwOut[r] = new BufferedWriter(new FileWriter(outputFile));
			r++;
		}
		fileName = path + outPrefix + "_gene.mcmc";
		outputFile = new File(fileName);
		bwOut[r] = new BufferedWriter(new FileWriter(outputFile));

		writeFileHeaderInformation(numberOfDomains);

	}

	// Method for writing header information on each output file
	public void writeFileHeaderInformation(int numberOfDomains) throws IOException {

		// Domain header
		String headerStr = String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", "iteration", "birthRate", "deathRate",
				"edgeMean", "edgeCV", "DLRModelDensity", "dataLikelihood", "overallDensity", "domainTree");
		for (int r = 0; r < numberOfDomains; r++) {
			bwOut[r].write(headerStr);
			bwOut[r].newLine();
			bwOut[r].flush();
		}

		// Gene Header
		headerStr = String.format("%s\t%s\t%s\t%s\t%s\t%s", "iteration", "birthRate", "deathRate", "DLModelDensity",
				"overallDensity", "geneTree");
		bwOut[numberOfDomains].write(headerStr);
		bwOut[numberOfDomains].newLine();

	}

	// This method will return domain level parameters and likelihood of models as
	// string record
	public String getDomainParametersStr(long iteration, Parameters pi, int domainIdx, StateLikelihood statelikelihood,
			boolean withLength, boolean withLeafNames) {

		LogDouble overall = new LogDouble(1.0);
		LogDouble domainDLRModel = statelikelihood.domainDLRModel;
		LogDouble domainSeqEvoModel = statelikelihood.domainSeqEvoModel[domainIdx];

		overall.mult(domainDLRModel);
		overall.mult(domainSeqEvoModel);

		StringBuilder lenStr = null;
		String rateStr, likelihoodStr, treeStr;
		StringBuilder sb = new StringBuilder();

		rateStr = String.format("%d\t%f\t%f\t%f\t%f", iteration, pi.domBDRates[domainIdx][0],
				pi.domBDRates[domainIdx][1], pi.domEdgeRates[domainIdx][0], pi.domEdgeRates[domainIdx][1]);

		likelihoodStr = String.format("%f\t%f\t%f ", domainDLRModel.getLogValue(), domainSeqEvoModel.getLogValue(),
				overall.getLogValue());

		if (withLength)
			treeStr = NewickWriter.getNewickString(pi.domainTree[domainIdx], pi.domainTree[domainIdx].bl,
					withLeafNames);
		else {
			treeStr = NewickWriter.getNewickString(pi.domainTree[domainIdx], null, withLeafNames);
			lenStr = new StringBuilder();
			for (int i = 0; i < pi.domainTree[domainIdx].nnodes; i++)
				lenStr.append(pi.domainTree[domainIdx].bl[i] + "\t");
			// lenStr.append(String.format("%13.8f \t", pi.domainTree[domainIdx].bl[i]));
		}

		sb.append(rateStr);
		sb.append("\t");
		sb.append(likelihoodStr);
		sb.append("\t");
		sb.append(treeStr);
		if (lenStr != null) {
			sb.append("\t");
			sb.append(lenStr);
		}

		return sb.toString();
	}

	// This method will return gene level parameters and likelihood of models as
	// string record
	public String getGeneParametersStr(long iteration, Parameters pi, StateLikelihood statelikelihood) {

		String rateStr, likelihoodStr, treeStr;
		StringBuilder sb = new StringBuilder();

		rateStr = String.format("%d\t%f\t%f", iteration, pi.geneBDRates[0], pi.geneBDRates[1]);

		likelihoodStr = String.format("%f\t%f", statelikelihood.geneDLModel.getLogValue(),
				statelikelihood.unnormalizedDensity.getLogValue());

		double branchTime = 0.0;
		double[] branchTimes = new double[pi.geneTree.nnodes];
		for (Node node : pi.geneTree.nodes) {
			if (!node.isRoot())
				branchTime = pi.geneTree.vt[node.parent.id] - pi.geneTree.vt[node.id];
			else
				branchTime = pi.geneTree.peakTime - pi.geneTree.vt[node.id];

			branchTimes[node.id] = branchTime * pi.speciesTree.leaf2TopTime;
		}

		treeStr = NewickWriter.getNewickString(pi.geneTree, branchTimes, true);

		sb.append(rateStr);
		sb.append("\t");
		sb.append(likelihoodStr);
		sb.append("\t");
		sb.append(treeStr);

		return sb.toString();
	}

	// This method will write domain level parameters and likelihood of models
	public void writeDomainParameters(long iteration, Parameters pi, int domainIdx, StateLikelihood statelikelihood,
			boolean withLength, boolean withLeafNames) throws IOException {
		String domainParamterStr = getDomainParametersStr(iteration, pi, domainIdx, statelikelihood, withLength,
				withLeafNames);
		bwOut[domainIdx].write(domainParamterStr);
		bwOut[domainIdx].newLine();
		bwOut[domainIdx].flush();
	}

	// This method will write gene level parameters and likelihood of models
	public void writeGeneParameters(long iteration, Parameters pi, StateLikelihood statelikelihood) throws IOException {
		String geneParamStr = getGeneParametersStr(iteration, pi, statelikelihood);
		bwOut[pi.numberOfDomains].write(geneParamStr);
		bwOut[pi.numberOfDomains].newLine();
		bwOut[pi.numberOfDomains].flush();
	}

	// This method will close output handler object
	public void close() throws IOException {
		for (int i = 0; i < bwOut.length; i++)
			bwOut[i].close();
	}

}

// treeStr = NewickWriter.getNewickString(pi.geneTree, null,true);
//////////
// private BufferedWriter [] bwTreeOut;
/*
 * 
 * //bwTreeOut[i].close();
 * 
 * // For only tree topology part of output bwTreeOut = new
 * BufferedWriter[numberOfDomains+1]; r=0; while(r<numberOfDomains) { fileName =
 * path + "top_dom_" + r + ".out"; outputFile = new File(fileName); bwTreeOut[r]
 * = new BufferedWriter(new FileWriter(outputFile)); r++; } fileName = path +
 * "top_gene.out"; outputFile = new File(fileName); bwTreeOut[r] = new
 * BufferedWriter(new FileWriter(outputFile));
 * 
 * 
 * 
 * // This method will write domain tree topology part only public void
 * writeDomainTopology(long iteration, Parameters pi,int domainIdx) throws
 * IOException {
 * 
 * StringBuilder sb = new StringBuilder();
 * 
 * sb.append( NewickWriter.getNewickTopologyString(pi.domainTree[domainIdx]) );
 * sb.append("\n");
 * 
 * bwTreeOut[domainIdx].write(sb.toString()); bwTreeOut[domainIdx].flush(); }
 * 
 * // This method will write gene tree topology part only public void
 * writeGeneTopology(long iteration, Parameters pi) throws IOException {
 * 
 * StringBuilder sb = new StringBuilder();
 * 
 * sb.append( NewickWriter.getNewickTopologyString(pi.geneTree) );
 * sb.append("\n");
 * 
 * bwTreeOut[pi.numberOfDomains].write(sb.toString());
 * bwTreeOut[pi.numberOfDomains].flush(); }
 * 
 * 
 */

/*
 * // Method for writing header information on each output file public void
 * writeFileHeaderInformation(int numberOfDomains) throws IOException {
 * 
 * // Domain header String headerStr = String.
 * format("%-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-15s \t %-15s \t %-15s \t %-15s"
 * , "iteration", "birthRate", "deathRate",
 * "edgeMean","edgeCV","DLRModelDensity","dataLikelihood","overallDensity",
 * "domainTree"); for(int r=0 ; r<numberOfDomains ; r++) {
 * bwOut[r].write(headerStr); bwOut[r].newLine(); bwOut[r].flush(); }
 * 
 * // Gene Header headerStr =
 * String.format("%-10s \t %-10s \t %-10s \t %-15s \t %-15s \t %-15s",
 * "iteration", "birthRate", "deathRate", "DLModelDensity",
 * "overallDensity","geneTree"); bwOut[numberOfDomains].write(headerStr);
 * bwOut[numberOfDomains].newLine();
 * 
 * }
 */

/*
 * 
 * // This method will write gene level parameters and likelihood of models
 * public void writeGeneParameters(long iteration, Parameters pi,StateLikelihood
 * statelikelihood) throws IOException {
 * 
 * String rateStr,likelihoodStr,treeStr; StringBuilder sb = new StringBuilder();
 * 
 * rateStr = String.format("%-10d \t %5.5f \t %5.5f", iteration,
 * pi.geneBDRates[0], pi.geneBDRates[1] );
 * 
 * likelihoodStr =
 * String.format("%-18.7f \t %-18.7f",statelikelihood.geneDLModel.getLogValue(),
 * statelikelihood.overall.getLogValue());
 * 
 * treeStr = NewickWriter.getNewickString(pi.geneTree, null,true);
 * 
 * sb.append(rateStr); sb.append("\t"); sb.append(likelihoodStr);
 * sb.append("\t"); sb.append(treeStr); sb.append("\n");
 * 
 * bwOut[pi.numberOfDomains].write(sb.toString());
 * bwOut[pi.numberOfDomains].flush(); }
 */
/*
 * // This method will write domain level parameters and likelihood of models
 * public void writeDomainParameters(long iteration, Parameters pi,int
 * domainIdx,StateLikelihood statelikelihood,boolean withLength,boolean
 * withLeafNames) throws IOException {
 * 
 * LogDouble overall = new LogDouble (1.0); LogDouble domainDLRModel =
 * statelikelihood.domainDLRModel; LogDouble domainSeqEvoModel =
 * statelikelihood.domainSeqEvoModel[domainIdx];
 * 
 * overall.mult(domainDLRModel); overall.mult(domainSeqEvoModel);
 * 
 * StringBuilder lenStr = null; String rateStr,likelihoodStr,treeStr;
 * StringBuilder sb = new StringBuilder();
 * 
 * 
 * rateStr = String.format("%-10d \t %5.5f \t %5.5f \t %5.5f \t %5.5f",
 * iteration, pi.domBDRates[domainIdx][0] , pi.domBDRates[domainIdx][1],
 * pi.domEdgeRates[domainIdx][0], pi.domEdgeRates[domainIdx][1]);
 * 
 * likelihoodStr = String.format("%-18.7f \t %-18.7f \t %-18.7f ",
 * domainDLRModel.getLogValue(), domainSeqEvoModel.getLogValue(),
 * overall.getLogValue());
 * 
 * 
 * if(withLength) treeStr =
 * NewickWriter.getNewickString(pi.domainTree[domainIdx],
 * pi.domainTree[domainIdx].bl,withLeafNames); else { treeStr =
 * NewickWriter.getNewickString(pi.domainTree[domainIdx], null,withLeafNames);
 * lenStr = new StringBuilder(); for(int i=0 ; i <
 * pi.domainTree[domainIdx].nnodes ; i++)
 * lenStr.append(String.format("%-10.7f \t", pi.domainTree[domainIdx].bl[i])); }
 * 
 * 
 * 
 * 
 * sb.append(rateStr); sb.append("\t"); sb.append(likelihoodStr);
 * sb.append("\t"); sb.append(treeStr); if(lenStr != null) { sb.append("\t");
 * sb.append(lenStr);} sb.append("\n");
 * 
 * bwOut[domainIdx].write(sb.toString()); bwOut[domainIdx].flush(); }
 */
