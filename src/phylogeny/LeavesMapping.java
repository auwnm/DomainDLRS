package phylogeny;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashMap;

public class LeavesMapping {

	public static LinkedHashMap<Integer, Integer> getMap(Tree guest, Tree host, Mapping guest2host) {

		int guestID, hostID;
		String guestLeaf, hostLeaf;

		LinkedHashMap<Integer, Integer> leafMap = new LinkedHashMap<Integer, Integer>();

		for (int i = 0; i < guest.nnodes; i++) {
			if (guest.nodes.get(i).isLeaf()) {
				guestID = guest.nodes.get(i).id;
				guestLeaf = guest.nodes.get(i).name;
				if (guest2host.map.containsKey(guestLeaf)) {
					hostLeaf = guest2host.map.get(guest.nodes.get(i).name);
					hostID = -1;
					for (int j = 0; j < host.nnodes; j++)
						if (host.nodes.get(j).isLeaf() && hostLeaf.equals(host.nodes.get(j).name))
							hostID = host.nodes.get(j).id;
					if (hostID == -1) {
						System.out.println("Error: Host Tree leaves name does not match to the mapping file names");
						return null;
					} else
						leafMap.put(guestID, hostID);
				} else {
					System.out.println("Error: Guest Tree leaves name does not match to the mapping file names");
					return null;
				}
			}
		}

		return leafMap;
	}

	// This function will read the mapping file and
	// return the hash map from guest leaves to host leaves.
	public static HashMap<Integer, Integer> Reader(String fileName, String[] hleaves, String[] gleaves)
			throws IOException {

		String line, guestID, hostID;
		int guestIndex, hostIndex;

		String[] tokens;
		HashMap<Integer, Integer> leaveMap = new HashMap<Integer, Integer>();

		File mappingsFile = new File(fileName);
		if (!mappingsFile.exists()) {
			System.out.println("Error: Cannot read sequence file " + fileName);
			return null;
		}

		BufferedReader brMappings = new BufferedReader(new FileReader(mappingsFile));

		while ((line = brMappings.readLine()) != null) {

			if (line.isEmpty())
				break;

			tokens = line.split("\\s+");
			guestID = tokens[0];
			hostID = tokens[1];

			guestIndex = hostIndex = -1;
			for (int i = 0; i < gleaves.length; i++) {
				if (guestID.equals(gleaves[i])) {
					guestIndex = i;
					break;
				}
			}
			if (guestIndex == -1) {
				System.out.println("Error: Guest Tree leaves name does not match to mapping file names");
				brMappings.close();
				return null;
			}

			for (int i = 0; i < hleaves.length; i++) {
				if (hostID.equals(hleaves[i])) {
					hostIndex = i;
					break;
				}
			}
			if (hostIndex == -1) {
				System.out.println("Error: Host Tree leaves name does not match to mapping file names");
				brMappings.close();
				return null;
			}

			leaveMap.put(guestIndex, hostIndex);
		}

		brMappings.close();
		return leaveMap;
	}

	public static int[] getLeavesMap(Tree guest, Tree host, String mappingFile) throws IOException {

		// Get leaves of host and guest trees.
		String[] hLeaves = host.getLeafNames(true);
		String[] gLeaves = guest.getLeafNames(true);

		// Indices Map induced from mapping file
		// From guest tree node list index ---> host tree node list index
		HashMap<Integer, Integer> indicesMap = Reader(mappingFile, hLeaves, gLeaves);

		// Leave mapping indexed by guest tree node ids ---> host tree node list indices
		int[] map = new int[guest.nnodes];
		for (int i = 0; i < guest.nnodes; i++) {
			if (guest.nodes.get(i).isLeaf()) {
				map[guest.nodes.get(i).id] = indicesMap.get(i);
			}
		}

		return map;
	}

	public static void show(Tree guest, Tree host, int[] map) {
		System.out.println("Leave Mapping:");
		for (int i = 0; i < guest.nnodes; i++) {
			if (guest.nodes.get(i).isLeaf()) {
				String geneName = String.format("%-30s", guest.nodes.get(i).name);
				String speciesName = String.format("%30s", host.nodes.get(map[guest.nodes.get(i).id]).name);
				System.out.println(geneName + "-->" + speciesName);
			}
		}

	}

}
