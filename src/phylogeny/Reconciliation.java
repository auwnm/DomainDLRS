package phylogeny;

import java.util.ArrayList;
import java.util.LinkedHashMap;

public class Reconciliation {

	public static final int EVENT_LEAF = 1;
	public static final int EVENT_SPEC = 2;
	public static final int EVENT_DUPL = 3;

	public Tree guest;
	public Tree host;

	public int[] leafMap;
	public int[] nodeMap;
	public int[] eventMap;

	public Tree guestPrime;
	public int[] nodeMapPrime;
	public int[] eventMapPrime;
	public int noImpliedNodes;
	public int[][] host2guestPrime;

	public Reconciliation(Tree guest, Tree host, Mapping guest2host) {
		this.guest = guest;
		this.host = host;
		this.leafMap = getMap(guest, host, guest2host);
		this.nodeMap = reconcile(guest, host, this.leafMap);
		this.eventMap = labelEvents(guest, this.nodeMap);
		this.noImpliedNodes = addImpliedSpecNodes(this.guest.copy(), this.host);
		this.host2guestPrime = getReverseMapPrime(host, this.guestPrime);
	}

	private int[][] getReverseMapPrime(Tree host, Tree guest) {

		int[][] reverseMap = new int[host.nnodes][guest.nnodes];

		for (int i = 0; i < host.nnodes; i++) {
			Node hnode = host.nodes.get(i);
			for (int j = 0; j < guest.nnodes; j++) {
				Node gnode = guest.nodes.get(j);
				if (this.nodeMapPrime[gnode.id] == hnode.id)
					reverseMap[hnode.id][gnode.id] = 1;
				else
					reverseMap[hnode.id][gnode.id] = -1;
			}
		}
		return reverseMap;
	}

	private int[] getMap(Tree guest, Tree host, Mapping guest2host) {

		int guestID, hostID;
		String guestLeaf, hostLeaf;

		int[] leafMap = new int[guest.nnodes];

		ArrayList<Node> guestLeaves = new ArrayList<Node>();
		for (Node node : guest.nodes) {
			if (node.isLeaf())
				guestLeaves.add(node);
			else
				leafMap[node.id] = -1;
		}

		ArrayList<Node> hostLeaves = new ArrayList<Node>();
		for (Node node : host.nodes) {
			if (node.isLeaf())
				hostLeaves.add(node);
		}

		for (Node gleaf : guestLeaves) {
			guestID = gleaf.id;
			guestLeaf = gleaf.name;
			if (guest2host.map.containsKey(guestLeaf)) {
				hostLeaf = guest2host.map.get(guestLeaf);
				hostID = -1;
				for (Node hleaf : hostLeaves) {
					if (hostLeaf.equals(hleaf.name)) {
						hostID = hleaf.id;
						break;
					}
				}
				if (hostID == -1) {
					System.out.println("Error: Host Tree leaves name does not match to the mapping file names");
					return null;
				} else
					leafMap[guestID] = hostID;
			} else {
				System.out.println("Error: Guest Tree leaves name does not match to the mapping file names");
				return null;
			}
		}
		return leafMap;
	}

	// Label events for each node in tree (assumes binary gene tree)
	private int[] labelEvents(Tree guest, int[] recon) {
		int[] eventMap = new int[guest.nnodes];
		for (Node node : guest.nodes) {
			if (node.nchildren == 0) {
				eventMap[node.id] = Reconciliation.EVENT_LEAF;
			} else {
				if (recon[node.id] == recon[node.children.get(0).id]
						|| recon[node.id] == recon[node.children.get(1).id])
					eventMap[node.id] = Reconciliation.EVENT_DUPL;
				else
					eventMap[node.id] = Reconciliation.EVENT_SPEC;
			}
		}
		return eventMap;
	}

	// Find Last Common Ancestor
	private Node treeLca(Tree htree, Node node1, Node node2) {
		int index1 = htree.nodeDepth[node1.id];
		int index2 = htree.nodeDepth[node2.id];

		while (index1 != index2) {
			if (index1 > index2) {
				node1 = node1.parent;
				index1 = htree.nodeDepth[node1.id];
			} else {
				node2 = node2.parent;
				index2 = htree.nodeDepth[node2.id];
			}
		}

		return node1;
	}

	// Recursively compute the reconciliation between guest and host tree (assumes
	// binary species tree)
	public void reconcile_recurse(Tree guest, Node node, Tree host, int[] recon) {
		// recurse
		for (int i = 0; i < node.nchildren; i++)
			reconcile_recurse(guest, node.children.get(i), host, recon);

		// post process
		if (node.nchildren > 0) {
			int sid1 = recon[node.children.get(0).id];
			int sid2 = recon[node.children.get(1).id];

			// this node's species is lca of children species
			Node lca = treeLca(host, host.nodes.get(sid1), host.nodes.get(sid2));
			recon[node.id] = lca.id;
		}
	}

	// TODO: implement more efficiently with post order traversal
	// reconcile a guest tree with a host tree
	private int[] reconcile(Tree guest, Tree host, int[] leafMap) {
		int[] recon = new int[guest.nnodes];

		// Label gene leaves with their species
		for (Node node : guest.nodes)
			if (node.isLeaf())
				recon[node.id] = leafMap[node.id];

		// Calling recursive part
		reconcile_recurse(guest, guest.root, host, recon);

		return recon;
	}

	// This method will add implied nodes in guest tree
	private int addImpliedSpecNodes(Tree gtree, Tree htree) {

		int addedNodes = 0;

		LinkedHashMap<Integer, Integer> nodeMapTemp = new LinkedHashMap<Integer, Integer>();
		LinkedHashMap<Integer, Integer> eventMapTemp = new LinkedHashMap<Integer, Integer>();

		for (Node node : gtree.nodes) {
			nodeMapTemp.put(node.id, this.nodeMap[node.id]);
			eventMapTemp.put(node.id, this.eventMap[node.id]);
		}

		// recurse
		int nnodes = gtree.nnodes;
		for (int i = 0; i < nnodes; i++) {
			Node node = gtree.nodes.get(i);
			// process this node and the branch above it

			// if no parent, then no implied speciation nodes above us
			// if (node->parent == NULL)
			// continue;

			// handle root node specially
			if (node.isRoot()) {
				// ensure root of gene tree properly reconciles to
				// root of species tree
				if (nodeMapTemp.get(node.id) == htree.root.id)
					continue;
				assert (gtree.root == node);
				// NOTE: root is not last node (may need to relax this condition)
				Node newnode = new Node();
				newnode.parent = null;
				newnode.isImplied = true;
				newnode.addChild(node);
				node.parent = newnode;

				gtree.root = gtree.addNode(newnode);
				nodeMapTemp.put(gtree.root.id, htree.root.id);
				eventMapTemp.put(gtree.root.id, Reconciliation.EVENT_SPEC);
				addedNodes++;

				/*
				 * // NOTE: root is not last node (may need to relax this condition) gtree.root
				 * = gtree.addNode(new Node()); gtree.root.children.add(node); node.parent =
				 * gtree.root; nodeMapTemp.put(gtree.root.id, htree.root.id);
				 * eventMapTemp.put(gtree.root.id, Reconciliation.EVENT_SPEC); addedNodes++;
				 */
			}

			// Determine starting and ending species
			Node sstart = htree.nodes.get(nodeMapTemp.get(node.id));
			Node send = htree.nodes.get(nodeMapTemp.get(node.parent.id));

			// The species path is too short to have implied speciations
			if (sstart == send)
				continue;

			Node parent = node.parent;

			// determine species path of this gene branch (node, node->parent)
			for (Node ptr = sstart.parent; ptr.id != send.id;) {
				// process ptr
				addSpecNode(node, ptr, gtree, nodeMapTemp, eventMapTemp);
				addedNodes++;
				node = node.parent;

				if (ptr.parent == null) {
					System.out.println(ptr.id);
					this.guest.show(false);

				}
				ptr = ptr.parent;

			}

			// Duplication Case:
			// Determine whether node.parent is a duplication
			// if so, send (a.k.a. species end) is part of species path
			if (eventMapTemp.get(parent.id) == Reconciliation.EVENT_DUPL) {
				addSpecNode(node, send, gtree, nodeMapTemp, eventMapTemp);
				addedNodes++;
			}
		}

		// Updating primes
		this.guestPrime = gtree;
		this.nodeMapPrime = new int[gtree.nnodes];
		this.eventMapPrime = new int[gtree.nnodes];

		for (Node node : gtree.nodes) {
			this.nodeMapPrime[node.id] = nodeMapTemp.get(node.id);
			this.eventMapPrime[node.id] = eventMapTemp.get(node.id);
		}
		return addedNodes;
	}

	// Insert new speciation node above node
	private void addSpecNode(Node node, Node snode, Tree gtree, LinkedHashMap<Integer, Integer> nodeMap,
			LinkedHashMap<Integer, Integer> eventMap) {
		Node newnode = new Node();
		Node parent = node.parent;

		// find index of node in parent's children
		int nodei = 0;
		for (; nodei < parent.nchildren; nodei++)
			if (parent.children.get(nodei).id == node.id)
				break;
		assert (nodei != parent.nchildren);

		// Insert new node into tree
		parent.children.set(nodei, newnode);
		newnode.parent = parent;
		newnode.addChild(node);
		node.parent = newnode;

		// Updating reconcliation and events info
		newnode.isImplied = true;
		gtree.addNode(newnode);
		nodeMap.put(newnode.id, snode.id);
		eventMap.put(newnode.id, Reconciliation.EVENT_SPEC);

	}

	public void show() {
		System.out.println();
		System.out.println("MPR Reconciliation:");
		int nnodes = guest.nnodes;
		for (int i = 0; i < nnodes; i++) {
			String guestID = String.format("%-2d", guest.nodes.get(i).id);
			String hostID = String.format("%-2d", nodeMap[guest.nodes.get(i).id]);
			int event = eventMap[guest.nodes.get(i).id];

			String showEvent = null;
			if (event == Reconciliation.EVENT_LEAF)
				showEvent = "Leaf";
			else if (event == Reconciliation.EVENT_SPEC)
				showEvent = "Speciation";
			else if (event == Reconciliation.EVENT_DUPL)
				showEvent = "Duplication";

			if (showEvent != null)
				System.out.println(guestID + " --> " + hostID + "\t" + showEvent);
			else
				System.out.println(guestID + " --> " + hostID);
		}
	}

	// This function will return reconciliation string for each node in guest tree
	public static String[] getReconStr(Tree guest, Tree host, int[] R) {

		String[] strArray = new String[guest.nnodes];

		for (int i = 0; i < guest.nnodes; i++) {
			StringBuilder sb = new StringBuilder();
			sb.append(String.format("%-2d", guest.nodes.get(i).id));
			sb.append(" --> ");
			sb.append(String.format("%-2d", host.nodes.get(R[guest.nodes.get(i).id]).id));
			strArray[i] = sb.toString();
		}

		return strArray;
	}

}
