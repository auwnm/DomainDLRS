package phylogeny;

import java.util.ArrayList;
import java.util.Arrays;

public class Tree {

	public String name; // Name of the tree
	public Node root; // Root of the tree (NULL if no nodes)
	public int nnodes; // Number of nodes in tree
	public ArrayList<Node> nodes; // Array of nodes (size = nnodes)
	public int[] nodeDepth; // Depth of each node indexed by node id
	public ArrayList<Node> postOrder; // List of tree nodes with post order arrangment;

	public double peakTime; // Peak point at stem edge
	public double leaf2TopTime; // Original peak time in my units
	public double[] vt; // Vertex times with respect to node id
	public double[] bl; // Branch lengths with respect to node id
	public int[] parray; // Parent of each node

	public double[] lengthCache; // hold the copy of branch lengths

	// Default constructor
	public Tree() {
		name = null;
		nnodes = 0;
		root = null;
		nodes = new ArrayList<Node>();
		nodeDepth = null;
		vt = null;
		bl = null;
		lengthCache = null;
		parray = null;
	}

	public Tree copy() {
		// Creating new tree object with the same name
		Tree tree = new Tree();
		tree.setName(this.name);

		// Allocate children with same attributes
		for (int i = 0; i < nnodes; i++) {
			Node node = new Node();
			node.id = this.nodes.get(i).id;
			node.name = this.nodes.get(i).name;
			node.originalLen = this.nodes.get(i).originalLen;
			node.isImplied = this.nodes.get(i).isImplied;
			tree.addNode(node);
		}

		// Copying time information
		if (this.vt != null)
			tree.vt = Arrays.copyOf(this.vt, this.vt.length);
		tree.peakTime = this.peakTime;

		// Store parent and child pointers
		int[] ptree = this.tree2ptree();
		for (int i = 0; i < nnodes; i++) {
			int parent = ptree[i];
			if (parent != -1) {
				Node parentnode = tree.nodes.get(parent);
				parentnode.addChild(tree.nodes.get(i));
				tree.nodes.get(i).parent = parentnode;
			} else {
				tree.nodes.get(i).parent = null;
			}
		}

		// Set root & depth of the tree
		tree.root = tree.nodes.get(nnodes - 1);
		if (this.nodeDepth != null)
			tree.setDepths(tree.root, 0);

		// Complete update
		if (this.parray != null)
			tree.parray = tree2ptree();

		if (this.postOrder != null) {
			tree.postOrder = new ArrayList<Node>();
			getTreePostOrder(this.postOrder, this.root);
		}

		// Validating tree structure and equality of two tree
		assert (tree.assertTree());
		assert (this.nnodes == tree.nnodes);
		assert (this.root.id == tree.root.id);
		for (int i = 0; i < this.nnodes; i++) {
			Node node1 = this.nodes.get(i);
			Node node2 = tree.nodes.get(i);

			assert (node1.id == node2.id);

			assert (this.vt[node1.id] == tree.vt[node2.id]);

			assert (node1.nchildren == node2.nchildren);

			if (node1.parent == null)
				assert (node2.parent == null);
			else
				assert (node1.parent.id == node2.parent.id);

			if (node1.isLeaf())
				assert (node2.isLeaf());
			else
				for (int j = 0; j < node1.nchildren; j++)
					assert (node1.children.get(j).id == node2.children.get(j).id);
		}
		return tree;
	}

	// Cache and restore methods to handle the length parameter
	public void cacheLengths() {
		if (lengthCache == null)
			this.lengthCache = Arrays.copyOf(this.bl, this.bl.length);
	}

	public void restoreLengths() {
		if (lengthCache != null)
			this.bl = lengthCache;
		this.lengthCache = null;
	}

	public void clearLengthsCache() {
		this.lengthCache = null;
	}

	// This will update the tree after perturbation (NNI, SPR)
	public void updateTree() {
		this.parray = tree2ptree();
		this.setDepths(this.root, 0);
		this.postOrder = new ArrayList<Node>();
		getTreePostOrder(this.postOrder, this.root);
	}

	// root tree by a new branch/node
	public void reroot(Node newroot, boolean onBranch) {
		// handle trivial case, newroot is root
		if (root == newroot || (onBranch && root.nchildren == 2
				&& (root.children.get(0) == newroot || root.children.get(1) == newroot)))
			return;

		// determine where to stop ascending
		Node oldroot = root;
		Node stop1 = null, stop2 = null;

		if (isRooted()) {
			stop1 = root.children.get(0);
			stop2 = root.children.get(1);
		} else {
			stop1 = root;
		}

		// start the reversal
		Node ptr1 = null, ptr2 = null;
		double nextDist = 0;
		double rootdist = 0.0;

		if (onBranch) {
			if (isRooted()) {
				// just need to stick current root somewhere else
				Node other = newroot.parent;
				rootdist = this.bl[stop1.id] + this.bl[stop2.id];

				oldroot.children.set(0, newroot);
				oldroot.children.set(1, other);
				newroot.parent = oldroot;
				// this.bl[newroot.id] /= 2.0;

				ptr1 = other;

				int oldchild = findval(ptr1.children, newroot);
				assert (oldchild != -1);

				// prepare for reversing loop
				ptr1.children.set(oldchild, oldroot);
				ptr2 = oldroot;
				nextDist = this.bl[newroot.id];
			} else {
				// need to add a new node to be root
				// TODO: not implemented
				// assert(0);
			}
		} else {
			if (isRooted()) {
				// need to remove the root node, and make tribranch
				// TODO: not implemented
				// assert(0);
			} else {
				// just need to swap node positions
				// TODO: not implemented
				// assert(0);
			}
		}

		// reverse parent child relationships
		while (ptr1 != stop1 && ptr1 != stop2) {
			int oldchild = findval(ptr1.children, ptr2);
			assert (oldchild != -1);

			Node next = ptr1.parent;

			// ptr1 is now fixed
			ptr1.children.set(oldchild, next);
			ptr1.parent = ptr2;

			// swap distances
			double tmpdist = this.bl[ptr1.id];
			this.bl[ptr1.id] = nextDist;
			nextDist = tmpdist;

			// move pointers
			ptr2 = ptr1;
			ptr1 = next;
		}

		// handle last two nodes
		if (stop2 != null) {
			// make stop1 parent of stop2
			if (stop2 == ptr1) {
				Node tmp = stop1;
				stop1 = ptr1;
				stop2 = tmp;
			}
			assert (ptr1 == stop1);

			int oldchild = findval(stop1.children, ptr2);
			stop1.children.set(oldchild, stop2);
			stop1.parent = ptr2;
			this.bl[stop1.id] = nextDist;
			stop2.parent = stop1;
			this.bl[stop2.id] = rootdist;
		} else {
			// assert(0);
		}

		// renumber nodes
		// - all leaves don't change numbers
		assert (root.id == (this.nnodes - 1));
	}

	int findval(ArrayList<Node> array, Node val) {
		for (int i = 0; i < array.size(); i++)
			if (array.get(i) == val)
				return i;
		return -1;
	}

	// This method will give arc time associated with a node
	public double getArcTime(Node node) {
		double time = 0.0;

		if (node.isRoot())
			time = getPeakTime() - this.vt[node.id];
		else
			time = this.vt[node.parent.id] - this.vt[node.id];

		return time;
	}

	// Sets the leaf names of the tree
	public void setLeafNames(String[] names, boolean leavesOnly) {
		for (int i = 0; i < nnodes; i++) {
			if (leavesOnly && !nodes.get(i).isLeaf())
				nodes.get(i).name = "";
			else
				nodes.get(i).name = names[i];
		}
	}

	// Gets leaf names of the nodes of a tree
	// Internal nodes are often named "" (empty string)
	public String[] getLeafNames(boolean leavesOnly) {
		String[] names = new String[nnodes];
		for (int i = 0; i < nnodes; i++)
			if (!leavesOnly || nodes.get(i).isLeaf())
				names[i] = nodes.get(i).name;
		return names;
	}

	// Gets number of leafs of this tree
	public int getNumberOfLeaves() {
		int leaves = 0;
		for (int i = 0; i < nnodes; i++)
			if (nodes.get(i).isLeaf())
				leaves++;
		return leaves;
	}

	// This method will return leave name based on its id
	public String getNodeName(int nodeID) {
		String name = null;
		for (int i = 0; i < nnodes; i++) {
			if (nodes.get(i).id == nodeID) {
				if (nodes.get(i).isLeaf())
					name = nodes.get(i).name;
				else
					name = String.format("%d", nodes.get(i).id);
			}
		}
		if (name == null)
			System.out.println("Error: could not be able to find the node of give id");

		return name;
	}

	// Gets names of the nodes of a tree
	// This differs from getLeafNames in that internal nodes will be named after
	// their name id (int) converted to a string.
	public String[] getNames() {
		String names[] = new String[nnodes];
		for (int i = 0; i < nnodes; i++) {
			if (nodes.get(i).isLeaf())
				names[i] = nodes.get(i).name;
			else
				names[i] = String.format("%d", nodes.get(i).id);
		}
		return names;
	}

	// Returns whether tree is rooted
	public boolean isRooted() {
		return (root != null && root.nchildren == 2);
	}

	// Returns the node with given id
	public Node getNode(int id) {
		for (int i = 0; i < this.nnodes; i++)
			if (this.nodes.get(i).id == id)
				return this.nodes.get(i);

		return null;
	}

	// Adds a node 'node' to the tree
	// This will also set the node's name id
	public Node addNode(Node node) {
		nodes.add(node);
		node.id = nodes.size() - 1;
		nnodes = nodes.size();
		return node;
	}

	public void getTreePostOrder(ArrayList<Node> postOrder, Node node) {
		if (node == null)
			node = this.root;

		// Recurse
		for (int i = 0; i < node.nchildren; i++)
			getTreePostOrder(postOrder, node.children.get(i));

		// Record post-process
		postOrder.add(node);

	}

	public ArrayList<Node> getTreePostList() {
		if (this.postOrder == null) {
			this.postOrder = new ArrayList<Node>();
			getTreePostOrder(this.postOrder, this.root);
		}
		return this.postOrder;
	}

	public void getTreePreOrder(ArrayList<Node> preOrder, Node node) {
		if (node == null)
			node = this.root;

		// Record pre-process
		preOrder.add(node);

		// Recurse
		for (int i = 0; i < node.nchildren; i++)
			getTreePostOrder(preOrder, node.children.get(i));

	}

	// create a parent tree from a tree object array
	public int[] tree2ptree() {
		int[] ptree = new int[this.nnodes];
		for (int i = 0; i < nnodes; i++) {
			if (nodes.get(i).parent != null)
				ptree[i] = nodes.get(i).parent.id;
			else
				ptree[i] = -1;
		}
		return ptree;
	}

	// This method will set depth of a tree
	public int setDepths(Node node, int index) {
		if (node == null)
			node = this.root;

		if (nodeDepth == null)
			nodeDepth = new int[this.nnodes];

		nodeDepth[node.id] = index;

		for (int i = 0; i < node.nchildren; i++)
			index = setDepths(node.children.get(i), ++index);

		return index;
	}

	// This method will set branch Length from original lengths
	public void setBranchLengthsFromOriginalLengths() {

		// Initialise branch Length arrays
		if (this.bl == null)
			this.bl = new double[this.nnodes];

		// Traverse each node
		for (int i = 0; i < this.nnodes; i++) {
			int idx = this.nodes.get(i).id;
			this.bl[idx] = this.nodes.get(i).originalLen;
		}

	}

	// This method will set times from original lengths
	public void setTimesFromOriginalLengths() {

		// Initialise the vertex time array
		if (vt == null)
			vt = new double[this.nnodes];

		// Setting vertex wise arc times
		double[] arcTimes = new double[this.nnodes];
		for (int i = 0; i < this.nnodes; i++) {
			int idx = this.nodes.get(i).id;
			arcTimes[idx] = this.nodes.get(i).originalLen;
			if (this.nodes.get(i).isLeaf())
				vt[idx] = 0.0;
			else
				vt[idx] = -1.0;
		}

		// Post Order Traversal for setting vertex times
		ArrayList<Node> postOrder = new ArrayList<Node>();
		getTreePostOrder(postOrder, this.root);
		for (Node node : postOrder) {
			if (!node.isRoot() && vt[node.id] != -1) {
				vt[node.parent.id] = vt[node.id] + arcTimes[node.id];
			}
		}

		// Setting peak time
		this.peakTime = vt[this.root.id] + arcTimes[this.root.id];
		this.leaf2TopTime = this.peakTime;

	}

	// Note: It will not touch arc time above root node
	public double[] getArcTimesFromVertexTimes() {

		double[] arctimes = new double[this.nnodes];

		ArrayList<Node> postOrder = new ArrayList<Node>();
		getTreePostOrder(postOrder, this.root);
		for (Node node : postOrder) {
			if (!node.isRoot() && vt[node.id] != -1)
				arctimes[node.id] = vt[node.parent.id] - vt[node.id];
		}

		return arctimes;
	}

	// This method will return the peak time of the tree
	public double getPeakTime() {
		return peakTime;
	}

	// This method will return the total arc time of the tree including root arc.
	public double getTotalTime() {
		if (this.vt == null)
			return Double.NaN;

		double totime = 0.0;
		for (Node node : this.nodes) {
			if (node.isRoot())
				totime += this.getPeakTime();
			else
				totime += (this.vt[node.parent.id] - this.vt[node.id]);
		}

		return totime;
	}

	// This method will return the normalized times (w.r.t root vertex time)
	public double[] getNormalizedTimes() {

		double[] normVT = new double[this.nnodes];

		double peakTime = getPeakTime();
		assert (peakTime > 0.0);

		// Rescale tree times so that peak node has vertex time 1.0
		for (int i = 0; i < this.nnodes; i++)
			normVT[i] = this.vt[i] / peakTime;

		return normVT;
	}

	// This method will return the normalized times (w.r.t root vertex time)
	public void setNormalizedTimes() {
		assert (this.peakTime > 0.0);

		// Rescale tree times so that peak has vertex time 1.0
		for (int i = 0; i < this.nnodes; i++)
			this.vt[i] /= this.peakTime;

		this.peakTime = 1.0;
	}

	public static Tree ptree2tree(int nnodes, int[] ptree) {
		Tree tree = new Tree();

		// allocate children
		for (int i = 0; i < nnodes; i++) {
			Node node = new Node();
			node.id = i;
			tree.addNode(node);
		}

		// store parent and child pointers
		for (int i = 0; i < nnodes; i++) {
			int parent = ptree[i];
			if (parent != -1) {
				Node parentnode = tree.nodes.get(parent);
				parentnode.addChild(tree.nodes.get(i));
				tree.nodes.get(i).parent = parentnode;
			} else {
				tree.nodes.get(i).parent = null;
			}
		}

		// set root
		tree.root = tree.nodes.get(nnodes - 1);
		tree.parray = ptree;
		assert (tree.assertTree());
		return tree;

	}

	// assert that the tree datastructure is self-consistent
	public boolean assertTree() {

		if (root == null) {
			System.out.println("root == null\n");
			return false;
		}
		if (nnodes != nodes.size()) {
			System.out.println("nnodes != nodes.size()\n");
			return false;
		}
		if (root.parent != null) {
			System.out.println("root.parent != null\n");
			return false;
		}
		/*
		 * if (root->name != nnodes - 1) { fprintf(stderr,
		 * "root->name != nnodes - 1\n"); return false; }
		 */

		boolean leaves = true;
		for (int i = 0; i < nnodes; i++) {
			// printf("assert %d\n", i);
			if (nodes.get(i) == null) {
				System.out.println("nodes[i] == null\n");
				return false;
			}

			// ids are correct
			if (nodes.get(i).id != i) {
				System.out.println("nodes[i].id != i\n");
				return false;
			}

			// do leaves come first
			if (nodes.get(i).isLeaf()) {
				if (!leaves) {
					System.out.println("!leaves\n");
					return false;
				}
			} else
				leaves = false;

			// check parent child pointers
			for (int j = 0; j < nodes.get(i).nchildren; j++) {
				// printf("assert %d %d\n", i, j);
				if (nodes.get(i).children.get(j) == null) {
					System.out.println("nodes[i]->children[j] == null\n");
					return false;
				}
				// printf("assert %d %d parent\n", i, j);
				if (nodes.get(i).children.get(j).parent != nodes.get(i)) {
					System.out.println("nodes[i]->children[j]->parent != nodes[i]\n");
					return false;
				}
			}
		}
		return true;
	}

	// This method will show nodes for debuging purposes
	public void show(boolean postOrder) {

		String child1 = new String();
		String child2 = new String();
		System.out.println("---------------------------");

		ArrayList<Node> nodeList;
		nodeList = this.nodes;

		if (postOrder) {
			nodeList = new ArrayList<Node>();
			getTreePostOrder(nodeList, this.root);
		}

		for (Node node : nodeList) {
			if (node.isLeaf())
				System.out.print("Leaf " + String.format("%3d", node.id) + " " + String.format("%-10s", node.name));
			else {
				System.out.print("Node " + String.format("%3d", node.id) + " ");
				child1 = child2 = null;
				for (int j = 0; j < node.nchildren; j++) {
					if (j != node.nchildren - 1)
						child1 = String.format("%-3d", node.children.get(j).id);
					else
						child2 = String.format("%-3d", node.children.get(j).id);
				}

				if (node.name != null)
					System.out.print(String.format("%-10s", node.name));

				String interName = "(" + child1 + "," + child2 + ")";
				System.out.print(String.format("%-10s", interName));
			}

			// if(this.at != null )
			// System.out.print(" at = " + String.format("%3.8f",this.at[node.id] ));
			// if(this.vt != null )
			// System.out.print(" vt = " + String.format("%3.8f",this.vt[node.id] ));

			if (node.isRoot())
				System.out.print("\t root node");
			System.out.println();
		}

	}

	// Setget the branch lengths arc and vertex times of the tree ordered w.r.t.
	// node ids
	public void setLengths(double[] lengths) {
		if (this.bl == null)
			this.bl = new double[nnodes];

		for (int i = 0; i < nnodes; i++)
			this.bl[i] = lengths[i];
	}

	public double[] getLengths() {
		return this.bl;
	}

	public void setVertexTimes(double[] vt) {
		this.vt = vt;
	}

	public double[] getVertexTimes() {
		return this.vt;
	}

	// Set name of the tree used for logging
	public void setName(String name) {
		this.name = name;
	}

	// Get the name of the tree
	public String getName() {
		return this.name;
	}

}

/*
 * 
 * //Note: It will not touch vertex time above root node
 * 
 * public boolean isUltrametric() {
 * 
 * boolean flag = true;
 * 
 * assert( vt != null ); //&& at != null
 * 
 * ArrayList<Node> postOrder = new ArrayList<Node>();
 * getTreePostOrder(postOrder,this.root); for(Node node: postOrder) { if(
 * !node.isRoot() && vt[node.id] != -1 ) { double differance =
 * vt[node.parent.id] - (vt[node.id] + at[node.id]); if(!(differance < 1e-10))
 * flag = false; } } return flag; }
 * 
 * 
 * 
 * // Note: It will not touch vertex time of (stem vertex) above root node
 * public void setVertexTimesFromArcTimes() {
 * 
 * if(vt == null) vt = new double [this.nnodes];
 * 
 * ArrayList<Node> postOrder = new ArrayList<Node>();
 * getTreePostOrder(postOrder,this.root);
 * 
 * for(int i=0;i<this.nnodes;i++) { int idx = this.nodes.get(i).id;
 * if(this.nodes.get(i).isLeaf()) vt[idx] = 0.0; else vt[idx] = -1; }
 * 
 * for(Node node: postOrder) { if( !node.isRoot() && vt[node.id] != -1 )
 * vt[node.parent.id] = vt[node.id] + at[node.id]; }
 * 
 * }
 * 
 * 
 * public double getLargerArcTime() {
 * 
 * if(this.at == null) return Double.NaN;
 * 
 * double max = Double.NEGATIVE_INFINITY; for(int i=0 ; i<this.nnodes ; i++)
 * if(this.at[i]>max) max = this.at[i];
 * 
 * return max; }
 * 
 */
