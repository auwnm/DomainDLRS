package phylogeny;

import java.util.ArrayList;

//A node in the phylogenetic tree
public class Node {

	public int id; // Node id (matches index in tree.nodes)
	public String name; // Node name (used mainly for leaves only)
	public Node parent; // Parent of this node
	public int nchildren; // Number of children
	public ArrayList<Node> children; // Array of child pointers (size = nchildren)
	public double originalLen; // Length measure associated with this node
	public boolean isImplied; // If this node is an implied host node

	// Default constructor
	public Node() {
		this.id = -1;
		this.originalLen = -1;
		this.parent = null;
		this.name = null;
		this.isImplied = false;
		this.nchildren = 0;
		this.children = new ArrayList<Node>();
	}

	public void addChild(Node node) {
		nchildren++;
		children.add(node);
		node.parent = this;
	}

	public boolean isLeaf() {
		return nchildren == 0;
	}

	public boolean isRoot() {
		return this.parent == null;
	}

}
