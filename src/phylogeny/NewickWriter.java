package phylogeny;

public class NewickWriter {

	public static boolean withLeafNames;

	// write out the newick notation of a tree
	public static void writeNewickNode(StringBuilder treeStr, Node node, double[] length, int depth) {

		String nodeName = withLeafNames ? node.name : String.format("v%d", node.id);
		if (node.nchildren == 0) {
			if (length != null)
				treeStr.append(String.format("%s:%f", nodeName, length[node.id]));
			else
				treeStr.append(String.format("%s", nodeName));
		} else {

			treeStr.append("(");
			for (int i = 0; i < node.nchildren - 1; i++) {
				writeNewickNode(treeStr, node.children.get(i), length, depth + 1);
				treeStr.append(",");
			}
			writeNewickNode(treeStr, node.children.get(node.nchildren - 1), length, depth + 1);

			if (nodeName != null)
				treeStr.append(")" + nodeName);
			else
				treeStr.append(")");

			if (depth >= 0 && length != null)
				treeStr.append(String.format(":%f", length[node.id]));
		}

	}

	/*
	 * // write out the newick notation of a tree public static void
	 * writeNewickNode(StringBuilder treeStr, Node node,double [] length, int depth)
	 * {
	 * 
	 * String nodeName = withLeafNames ? node.name : String.format("v%d",node.id) ;
	 * if (node.nchildren == 0) { if(length!=null) treeStr.append(
	 * String.format("%s:%10.7f", nodeName, length[node.id]) ); else
	 * treeStr.append(String.format("%s", nodeName)); } else {
	 * 
	 * treeStr.append("("); for (int i=0; i<node.nchildren - 1; i++) {
	 * writeNewickNode(treeStr,node.children.get(i), length, depth+1);
	 * treeStr.append(","); }
	 * writeNewickNode(treeStr,node.children.get(node.nchildren-1),length,depth+1);
	 * treeStr.append(")"+ nodeName ); if (depth > 0 && length != null )
	 * treeStr.append( String.format(":%f", length[node.id] ) ); }
	 * 
	 * 
	 * }
	 */
	// Output the newick notation of a tree in string
	public static String getNewickString(Tree tree, double[] length, boolean withLeafNames) {
		NewickWriter.withLeafNames = withLeafNames;
		StringBuilder treeStr = new StringBuilder();
		writeNewickNode(treeStr, tree.root, length, 0);
		treeStr.append(";");
		return treeStr.toString();
	}

	// Output only topology part of the tree
	public static String getNewickTopologyString(Tree tree) {
		NewickWriter.withLeafNames = true;
		double[] length = null;
		StringBuilder treeStr = new StringBuilder();
		writeNewickNode(treeStr, tree.root, length, 0);
		treeStr.append(";");
		return treeStr.toString();
	}

}
