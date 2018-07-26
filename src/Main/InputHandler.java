package Main;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedHashMap;

public class InputHandler {

	public static int domainMapCounter = 0;
	public static int domainMSACounter = 0;
	public static LinkedHashMap<String, Integer> parameterMap = new LinkedHashMap<String, Integer>();

	public static Input readInputArguments(String[] args) throws IOException {
		Input input = null;
		if (args.length < 1 || args[0].equalsIgnoreCase("-h") || args[0].equalsIgnoreCase("--help")) {
			System.out.println("================================================================================\n"
					+ "DomainDLRS is a MCMC based application for studying the domain evolution within \n"
					+ "gene & species trees developed at KTH Computational Biology Department with     \n"
					+ "Science for Life Laboratory (SciLifeLab) in Stockholm, Sweden.\n\n"
					+ "Releases, source code and tutorial: https://github.com/auwnm/DomainDLRS\n\n"
					+ "License: DomainDLRS is available under the New XYZ License.\n");
			System.out.println("Usage: java -jar domainDLRS.jar -p <ParameterFile> \n");
			System.out.println("Help : java -jar domainDLRS.jar -h                 \n");
			System.out.println("================================================================================\n");
		} else if (args.length > 0 && args[0].equals("-p")) {
			input = readParamtereFile(args[1]);
		} else
			System.out.println("Help : java -jar domainDLRS.jar -h\n");

		return input;
	}

	public static Input readParamtereFile(String fileName) throws IOException {

		File infile = new File(fileName);
		if (!infile.exists())
			System.out.println("Error: Parameter file does not exist with this name ...");
		BufferedReader reader = new BufferedReader(new FileReader(infile));

		String line;
		LoadParameterMap();
		Input input = new Input();

		while ((line = reader.readLine()) != null) {

			// Dealing with comments
			if (line.startsWith("#") || line.isEmpty())
				continue;
			String[] tokens = line.split("=");

			// Expected format
			if (tokens.length == 2) {
				String parameter = tokens[0].trim();
				String value = tokens[1].trim();
				if (parameterMap.containsKey(parameter)) {
					if (value.equalsIgnoreCase("null"))
						AssignParameterValue(parameterMap.get(parameter), null, input);
					else
						AssignParameterValue(parameterMap.get(parameter), value, input);
				}
			}

		}

		reader.close();
		return input;
	}

	// Method will assign the values
	private static void AssignParameterValue(int paramNo, String value, Input input) {

		switch (paramNo) {

		case 0:
			input.path = value;
			break;
		case 1:
			input.numberOfDomains = Integer.parseInt(value);
			input.domMappingFiles = new String[input.numberOfDomains];
			input.domAlignmentFiles = new String[input.numberOfDomains];
			break;
		case 2:
			input.speciesTreeFile = value;
			break;
		case 3:
			input.geneMapFile = value;
			break;
		case 4:
			input.geneTree = value;
			break;
		case 5:
			input.domMappingFiles[domainMapCounter++] = value;
			break;
		case 6:
			input.domAlignmentFiles[domainMSACounter++] = value;
			break;
		case 7:
			input.substitutionModel = value;
			break;
		case 8:
			TuningParameters.mcmcMaxIteration = Integer.parseInt(value);
			break;
		case 9:
			TuningParameters.thinningFactor = Integer.parseInt(value);
			break;
		case 10:
			TuningParameters.useNativeCode = Boolean.parseBoolean(value);
			break;
		case 11:
			if (value == null)
				input.outPrefix = "result";
			else
				input.outPrefix = new String(value);
			break;
		case 12:
			TuningParameters.fixedGeneTree = Boolean.parseBoolean(value);
			break;
		case 13:
			TuningParameters.runType = new String(value);
			break;
		default:

		}
	}

	// Method with Parameter specs
	private static void LoadParameterMap() {
		parameterMap.put("path", 0);
		parameterMap.put("numberOfDomains", 1);
		parameterMap.put("speciesTreeFile", 2);
		parameterMap.put("geneMap", 3);
		parameterMap.put("geneTree", 4);
		parameterMap.put("domainMap", 5);
		parameterMap.put("domainMSA", 6);
		parameterMap.put("substitutionModel", 7);
		parameterMap.put("maxIteration", 8);
		parameterMap.put("thinningFactor", 9);
		parameterMap.put("useNativeCode", 10);
		parameterMap.put("outPrefix", 11);
		parameterMap.put("fixedGeneTree", 12);
		parameterMap.put("runType", 13);
	}

}
