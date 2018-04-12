package seqEvolution;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class FastaRW {

	public static Sequences readFastaAlign(String fileName) throws IOException {

		File infile = new File(fileName);
		if (!infile.exists())
			return null;

		BufferedReader reader = new BufferedReader(new FileReader(infile));
		Sequences sequences = new Sequences();

		String line;
		boolean lastSequence = false;
		StringBuilder seqStr = new StringBuilder();
		StringBuilder seqID = new StringBuilder();

		while ((line = reader.readLine()) != null) {
			if (line.startsWith(">")) {
				if (lastSequence)
					sequences.append(seqID.toString(), seqStr.toString());
				seqStr = new StringBuilder();
				seqID = new StringBuilder();
				seqID.append(line.substring(1).trim());
				lastSequence = false;
			} else {
				seqStr.append(line.trim());
				lastSequence = true;
			}
		}
		if (lastSequence && seqStr.length() != 0)
			sequences.append(seqID.toString(), seqStr.toString());

		reader.close();
		return sequences;
	}

	public static boolean write(String filename, Sequences sequences) throws IOException {
		File outfile = new File(filename);
		BufferedWriter writer = new BufferedWriter(new FileWriter(outfile));
		for (String key : sequences.seqMap.keySet()) {
			writer.write(">" + key);
			writer.newLine();
			writer.write(sequences.seqMap.get(key));
			writer.newLine();
		}
		writer.flush();
		writer.close();
		return true;
	}
}
