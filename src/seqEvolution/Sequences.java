package seqEvolution;

import java.util.LinkedHashMap;

public class Sequences {

	public int nseqs;
	public int seqlen;

	// Use of LinkedHashMap,
	// It will iterate in the order in which the entries were put into the map.
	public LinkedHashMap<String, String> seqMap;

	// Default Constructor
	public Sequences() {
		nseqs = seqlen = -1;
		seqMap = new LinkedHashMap<String, String>();
	}

	// This method will append new sequence
	public void append(String seqIdentifier, String seqString) {
		assert (this.seqlen == seqString.length() || this.seqlen == -1);
		assert (!this.seqMap.containsKey(seqIdentifier));
		seqMap.put(seqIdentifier, seqString);
		seqlen = seqString.length();
		nseqs = seqMap.size();
	}

}
