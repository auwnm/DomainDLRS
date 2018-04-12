package phylogeny;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;

public class NewickReader {

	private static int pos = 0;
	private static int depth = 0;
	private static String treeStr = null;
	private static String token = null;

	static char readChar() {
		char chr;

		// Indicate End Of String (EOS)
		if (pos == treeStr.length())
			return '\0';

		chr = treeStr.charAt(pos);
		pos++;

		// keep track of parent depth
		if (chr == '(')
			depth++;
		if (chr == ')')
			depth--;

		return chr;
	}

	static char readUntil(String stops) {
		char chr;
		token = "";

		while (true) {
			chr = readChar();

			if (chr == '\0')
				return chr;

			// compare char to stop characters
			for (int i = 0; i < stops.length(); i++) {
				if (chr == stops.charAt(i))
					return chr;
			}
			token += chr;
		}
	}

	static Node readNewickNode(Tree tree, Node parent) {
		char char1;
		Node node;

		// Read first Character
		if ((char1 = readChar()) == '\0') {
			System.out.println("Error: Not a valid Newick tree string ...");
			return null;
		}

		if (char1 == '(') {

			// read internal node
			int depth2 = depth;
			node = tree.addNode(new Node());
			if (parent != null)
				parent.addChild(node);

			// read all child nodes at this depth
			while (depth == depth2) {
				Node child = readNewickNode(tree, node);
				if (child == null)
					return null;
			}

			// Read name
			char ch = readUntil(":),;");

			// In case internal nodes have names
			if (token.length() != 0) {
				node.name = new String(token);
			}

			// read distance for this node
			if (ch == ':') {
				ch = readUntil("),;");

				if (token.length() != 0)
					node.originalLen = Double.parseDouble(token);

				if (ch == '\0')
					return null;
			}

			return node;

		} else {

			// read leaf
			node = tree.addNode(new Node());
			if (parent != null)
				parent.addChild(node);

			char ch;

			// Read name
			ch = readUntil(":),");
			if (ch == '\0')
				return null;
			token = char1 + token;
			node.name = token;

			// Read distance for this node
			if (ch == ':') {
				ch = readUntil(",)");

				if (token.length() != 0)
					node.originalLen = Double.parseDouble(token);

				if (ch == '\0')
					return null;
			}

			return node;
		}
	}

	/*
	 * static Node readNewickNode(Tree tree, Node parent) { char char1; Node node;
	 * 
	 * //Read first Character if( (char1 = readChar()) == '\0' ) {
	 * System.out.println("Error: Not a valid Newick tree string ..."); return null;
	 * }
	 * 
	 * 
	 * if (char1 == '(') {
	 * 
	 * // read internal node int depth2 = depth; node = tree.addNode(new Node()); if
	 * (parent != null) parent.addChild(node);
	 * 
	 * // read all child nodes at this depth while (depth == depth2) { Node child =
	 * readNewickNode(tree, node); if (child==null) return null; }
	 * 
	 * // read distance for this node char ch = readUntil( "):," ); if(ch == ':') {
	 * if( readUntil( "):,;" ) == '\0' ) return null; node.originalLen =
	 * Double.parseDouble(token); }
	 * 
	 * return node;
	 * 
	 * } else {
	 * 
	 * // read leaf node = tree.addNode(new Node()); if (parent != null)
	 * parent.addChild(node);
	 * 
	 * // Read name if ( readUntil( ":)," ) == '\0' ) return null; token = char1 +
	 * token; node.name = token;
	 * 
	 * // read distance for this node if( readUntil( ":)," ) == '\0' ) return null;
	 * node.originalLen = Double.parseDouble(token);
	 * 
	 * return node; } }
	 * 
	 */

	// This function will read file and buffer its contents in string format
	// This will also remove all newlines,carriage retruns,spaces and tabs from file
	// contents.
	static String readFileAsString(String filePath) throws IOException {
		StringBuffer fileData = new StringBuffer();
		BufferedReader reader = new BufferedReader(new FileReader(filePath));
		char[] buf = new char[1024];
		int numRead = 0;
		while ((numRead = reader.read(buf)) != -1) {
			String readData = String.valueOf(buf, 0, numRead);
			fileData.append(readData);
		}
		reader.close();
		String fileContents = fileData.toString().replaceAll("\\n|\\r|\\s|\\t", "");
		return fileContents;
	}

	public static Tree readNewickTreeFile(String fileName) throws IOException {
		String NewickTreeStr = NewickReader.readFileAsString(fileName);
		return readNewickTreeStr(NewickTreeStr);
	}

	/*
	 * Adapted from Rasmussen et. al. 2011 Parent Array Format (ptree)
	 * 
	 * 
	 * 4 / \ 3 \ / \ \ / \ \ 0 1 2
	 * 
	 * Then the parent tree array representation is ptree = [3, 3, 4, 4, -1] such
	 * that ptree[node's id] = node's parent's id
	 * 
	 * In addition, the following must be true 1. tree must be binary: n leaves, n-1
	 * internal nodes 2. leaves must be numbered 0 to n-1 3. internal nodes are
	 * numbered n to 2n-2 4. root must be numbered 2n-2 5. the parent of root is -1
	 * 6. the length of ptree is 2n-1
	 */

	// Main newick tree reader function
	public static Tree readNewickTreeStr(String NewickTreeStr) throws IOException {

		// Initializing the variables
		pos = 0;
		depth = 0;
		token = new String();
		treeStr = NewickTreeStr;
		Tree tree = new Tree();
		tree.root = readNewickNode(tree, null);

		// tree.show(false);

		// Changing root position
		tree.nodes.set(tree.root.id, tree.nodes.get(tree.nnodes - 1));
		tree.nodes.set(tree.nnodes - 1, tree.root);
		tree.nodes.get(tree.root.id).id = tree.root.id;
		tree.nodes.get(tree.nnodes - 1).id = tree.nnodes - 1;

		// tree.show(false);

		// Sorting the nodes list w.r.t. leaves and internal nodes
		Collections.sort(tree.nodes, new nodeComparator());

		// Assigning new order to nodes
		for (int i = 0; i < tree.nnodes; i++)
			tree.nodes.get(i).id = i;

		// For debuging purposes
		// At this stage result is not same as Rasmussen & Kellis
		// But their constraints has been met.
		// tree.show();

		return tree;
	}

}
