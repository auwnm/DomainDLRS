package posteriorAnalysis;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import phylogeny.NewickReader;
import phylogeny.NewickWriter;
import phylogeny.Node;
import phylogeny.Tree;

public class PosteriorReader {

	public static boolean mrBayesPosterior = false;

	// Reading MrBayes Posterior
	// Note onlyMapTree means the tree at the top of posterior made by MrBayes
	public static ArrayList<String> readMrBayesPosterior(String fileName, boolean onlyMapTree) throws IOException {

		ArrayList<String> nwkStrList = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new FileReader(fileName));

		String line;
		HashMap<Integer, String> map = new HashMap<Integer, String>();
		boolean translateTable = false;
		boolean treeSection = false;

		ArrayList<Double> posteriors = new ArrayList<Double>();
		ArrayList<String> treeStrList = new ArrayList<String>();

		while ((line = br.readLine()) != null) {

			String[] token = line.trim().split("\\s+");

			if (token.length == 0)
				continue;

			if (token[0].equalsIgnoreCase("translate")) {
				translateTable = true;
				continue;
			}

			if (token[0].equalsIgnoreCase("tree")) {
				translateTable = false;
				treeSection = true;
			}
			if (translateTable) {
				int id = Integer.parseInt(token[0]);
				String name = token[1].substring(0, token[1].length() - 1);
				map.put(id, name);
			}

			if (token[0].equalsIgnoreCase("end;"))
				break;

			if (treeSection) {

				String treeStr = token[11];
				double P = Double.parseDouble(token[7].split("]")[0]);

				// Top tree in the MrBayes posterior
				if (onlyMapTree) {
					treeStrList.add(treeStr);
					break;
				}

				// Taking MrBayes Posterior
				if (mrBayesPosterior) {
					double p = Double.parseDouble(token[4].split(",")[0]);
					posteriors.add(p);
					treeStrList.add(treeStr);
				} else if (P < 0.960) // Taking all trees fall in 95% of the posterior
					treeStrList.add(treeStr);

			}
		}
		br.close();

		// Preparing Posterior ...
		for (int i = 0; i < treeStrList.size(); i++) {
			// System.out.println(postTree);
			Tree tree = NewickReader.readNewickTreeStr(treeStrList.get(i));
			ArrayList<Node> nodes = tree.nodes;
			for (Node node : nodes) {
				if (node.isLeaf()) {
					String nodeName = node.name;
					node.name = map.get(Integer.parseInt(nodeName));
				}
			}
			// tree.show(false);
			// System.out.println);
			String treeStr = NewickWriter.getNewickTopologyString(tree);

			if (mrBayesPosterior) {
				double p = posteriors.get(i);
				if (p == 0)
					continue; // Note: can neglect cases with posterior=0
				int numberOfTimes = (int) (100.0 * p);

				if (numberOfTimes <= 1)
					nwkStrList.add(treeStr);
				else
					for (int jj = 1; jj < numberOfTimes; jj++)
						nwkStrList.add(treeStr);

			} else
				nwkStrList.add(treeStr);

		}

		return nwkStrList;
	}

	// Reading JPrime Posterior
	// Note onlyMapTree means the tree with "max posterior probability"
	public static ArrayList<String> readJPrimePosterior(String fileName, boolean onlyMapTree) throws IOException {

		int treeCol = 9;
		int posteriorCol = 2;
		boolean firstLine = true;

		String mapTree = null;
		double maxAposterior = Double.NEGATIVE_INFINITY;
		ArrayList<String> nwkStrList = new ArrayList<String>();

		String line;
		BufferedReader br = new BufferedReader(new FileReader(fileName));
		while ((line = br.readLine()) != null) {

			// Leaving First header Line
			if (firstLine) {
				firstLine = false;
				continue;
			}

			// Leaving empty lines
			if (line.length() == 0)
				break;

			String[] token = line.split("\\t");
			double posterior = Double.parseDouble(token[posteriorCol - 1]);
			String tree = token[treeCol - 1];

			if (posterior > maxAposterior) {
				mapTree = new String(tree);
				maxAposterior = posterior;
			}

			if (!onlyMapTree)
				nwkStrList.add(tree.trim());

		}
		br.close();

		if (onlyMapTree)
			nwkStrList.add(mapTree.trim());

		return nwkStrList;
	}

	// This method will read domainDLRS .inf file in specified format for hill
	// climbing solution
	public static String[] readHCFile(String fileName, int numberOfDomain) throws IOException {

		String[] data = new String[numberOfDomain + 1];
		BufferedReader br = new BufferedReader(new FileReader(fileName));

		int idx = 0;
		String line;
		boolean flag = true;
		while ((line = br.readLine()) != null) {

			String[] token = line.split("\\t");

			if (token.length == 0)
				continue;

			if (line.startsWith("# Best encountered state:")) {
				flag = false;
			}

			if (flag)
				continue;

			if (!token[token.length - 1].endsWith(";"))
				continue;
			else
				data[idx++] = token[token.length - 1];

		}
		br.close();

		return data;
	}

	// This method will read domainDLRS gene file in specified format.
	public static ArrayList<Sample> readGeneRunFile(String fileName, boolean header) throws IOException {
		ArrayList<Sample> samples = new ArrayList<Sample>();

		BufferedReader br = new BufferedReader(new FileReader(fileName));

		int i;
		String line;
		i = 0;
		while ((line = br.readLine()) != null) {

			// Leaving first line as header
			if (header && i == 0) {
				i++;
				continue;
			}

			// Empty line
			if (line.length() == 0)
				break;

			// Parsing line
			Sample sample = new Sample();
			String[] token = line.split("\\t");
			sample.iteration = Integer.parseInt(token[0]);
			sample.birthRate = Double.parseDouble(token[1]);
			sample.deathRate = Double.parseDouble(token[2]);
			sample.DLRModelDensity = Double.parseDouble(token[3]);
			sample.overallDensity = Double.parseDouble(token[4]);
			sample.tree = token[5].trim();

			// Adding sample
			samples.add(sample);
			i++;

		}
		br.close();

		return samples;
	}

	// This method will read domainDLRS domain file in a specified format.
	public static ArrayList<Sample> readDomainRunFile(String fileName, boolean header) throws IOException {

		ArrayList<Sample> samples = new ArrayList<Sample>();
		BufferedReader br = new BufferedReader(new FileReader(fileName));

		int i;
		String line;
		i = 0;
		while ((line = br.readLine()) != null) {

			// Leaving first line as header
			if (header && i == 0) {
				i++;
				continue;
			}

			// Empty line
			if (line.length() == 0)
				break;

			// Parsing line
			Sample sample = new Sample();
			String[] token = line.split("\\t");
			sample.iteration = Integer.parseInt(token[0]);
			sample.birthRate = Double.parseDouble(token[1]);
			sample.deathRate = Double.parseDouble(token[2]);
			sample.edgeMean = Double.parseDouble(token[3]);
			sample.edgeCV = Double.parseDouble(token[4]);
			sample.DLRModelDensity = Double.parseDouble(token[5]);
			sample.dataLikelihood = Double.parseDouble(token[6]);
			sample.overallDensity = Double.parseDouble(token[7]);
			sample.tree = token[8].trim();

			// Adding sample
			samples.add(sample);
			i++;

		}
		br.close();

		return samples;

	}

	// This method will read domainDLRS gene file in specified format.
	public static ArrayList<Sample> readSingleColFile(String fileName, boolean header) throws IOException {
		ArrayList<Sample> samples = new ArrayList<Sample>();

		BufferedReader br = new BufferedReader(new FileReader(fileName));

		int i;
		String line;
		i = 0;
		while ((line = br.readLine()) != null) {

			// Leaving first line as header
			if (header && i == 0) {
				i++;
				continue;
			}

			// Empty line
			if (line.length() == 0)
				break;

			// Parsing line
			Sample sample = new Sample();
			sample.iteration = i;
			sample.tree = line.trim();

			// Adding sample
			samples.add(sample);
			i++;

		}
		br.close();

		return samples;
	}

}
