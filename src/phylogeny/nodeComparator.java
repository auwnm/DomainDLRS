package phylogeny;

import java.util.Comparator;

public class nodeComparator implements Comparator<Node> {
	@Override
	public int compare(Node n1, Node n2) {
		if (n1.isLeaf()) {
			if (n2.isLeaf())
				return 0;
			else
				return -1;
		} else {
			if (n2.isLeaf())
				return 1;
			else {

				if (n1.id > n2.id)
					return 1;
				else
					return -1;

			}
			// return 0;
		}
	}
}

/*
 * 
 * @Override public int compare(Node n1, Node n2) { if ( n1.isLeaf() ) {
 * if(n2.isLeaf()) return 0; else return -1; } else{ if (n2.isLeaf()) return 1;
 * else return 0; } }
 */