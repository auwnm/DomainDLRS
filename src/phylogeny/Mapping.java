package phylogeny;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

public class Mapping {

	public LinkedHashMap<String, String> map;

	// Constructor
	public Mapping() {
		map = new LinkedHashMap<String, String>();
	}

	// Reading Mapping file
	public boolean readMapping(String fileName) throws IOException {

		String line, guestID, hostID;
		String[] tokens;
		File mappingsFile = new File(fileName);

		if (!mappingsFile.exists()) {
			map = null;
			return false;
		}

		BufferedReader brMappings = new BufferedReader(new FileReader(mappingsFile));
		while ((line = brMappings.readLine()) != null) {
			if (line.isEmpty())
				break;
			tokens = line.split("\\s+");
			guestID = tokens[0];
			hostID = tokens[1];
			map.put(guestID, hostID);
		}
		brMappings.close();
		return true;
	}

	// Show the mappings
	public void showMapping() {
		for (Map.Entry<String, String> entry : map.entrySet()) {
			String key = entry.getKey();
			String value = entry.getValue();
			System.out.println(key + "\t" + value);
		}
	}

}
