package phylogeny;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;

import org.forester.archaeopteryx.AptxUtil;
import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.archaeopteryx.Archaeopteryx;
import org.forester.archaeopteryx.Configuration;
import org.forester.archaeopteryx.Options;
import org.forester.archaeopteryx.TreeColorSet;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.BranchData;
import org.forester.phylogeny.data.BranchWidth;
import org.forester.phylogeny.data.NodeVisualization.NodeFill;
import org.forester.phylogeny.data.NodeVisualization.NodeShape;

import Main.mainSystem;

public class visTree {

	public static void main(String[] args) throws IOException {

		boolean putLeafName = false;

		String speciesTreeFile = mainSystem.CWD + "hostTree.nwk";
		String outfile = mainSystem.CWD + "Atree.png";

		Tree speciesTree = NewickReader.readNewickTreeFile(speciesTreeFile);
		speciesTree.setName("speciesTree");

		// Getting Phylogeny object of the tree.
		final Phylogeny phy = getPhylogeny(speciesTree, putLeafName, false);

		// Setting up a configuration object.
		final Configuration config = new Configuration();
		config.putDisplayColors(TreeColorSet.BACKGROUND, new Color(255, 255, 255));
		config.putDisplayColors(TreeColorSet.BRANCH, new Color(0, 0, 0));
		config.putDisplayColors(TreeColorSet.TAXONOMY, new Color(0, 0, 0));

		// Speciations = Red Duplications = Green Color
		config.putDisplayColors(TreeColorSet.SPECIATION, new Color(255, 0, 0));
		config.putDisplayColors(TreeColorSet.DUPLICATION, new Color(0, 255, 0));
		config.putDisplayColors(TreeColorSet.DUPLICATION_OR_SPECATION, new Color(0, 0, 255));

		// Node Shape
		config.setDefaultNodeShape(NodeShape.CIRCLE);
		config.setDefaultNodeFill(NodeFill.SOLID);
		config.setShowDefaultNodeShapesExternal(true);
		config.setShowDefaultNodeShapesInternal(true);
		config.setDefaultNodeShapeSize((short) 8);
		config.setBaseFontSize(14);

		// Branch related
		config.setColorizeBranches(true);
		config.setUseBranchesWidths(true);

		// Tree Style
		config.setPhylogenyGraphicsType(Options.PHYLOGENY_GRAPHICS_TYPE.ROUNDED);
		config.setDisplayNodeNames(true);
		config.setTaxonomyColorize(true);
		config.setDisplayTaxonomyCode(false);

		// Image Resolution
		int[] resolution = { 650, 650 };

		// Writing to a graphics file.
		AptxUtil.writePhylogenyToGraphicsFile(phy, new File(outfile), resolution[0], resolution[1],
				GraphicsExportType.PNG, config);

	}

	public static void showRecon(Tree guest, Tree host, String filePath) throws IOException {

		// Image names
		String img1 = "ImageHost.png";
		String img2 = "ImageGuest.png";

		// Image Resolution
		int[] resolution = { 650, 650 };

		// Drawing trees
		drawHostTree(host, filePath + img1, true, resolution);
		drawGuestTree(guest, filePath + img2, false, resolution);

		// Combine them into single image
		File path = new File(filePath);

		// Loading the source images
		BufferedImage image = rotateImage(ImageIO.read(new File(path, img1)), 90);
		BufferedImage overlay = rotateImage(ImageIO.read(new File(path, img2)), 90);

		// create the new image, canvas size is the max. of both image sizes
		int w = 1300;
		int h = 700;
		BufferedImage combined = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB);

		// Setting Background color
		Graphics g = combined.getGraphics();
		g.setColor(Color.WHITE);
		g.fillRect(0, 0, w, h);

		// paint both images, preserving the alpha channels
		g.drawImage(image, 0, 0, null);
		g.drawImage(overlay, 650, 0, null);

		// Deleting temporary files
		File temp1 = new File(filePath + img1);
		File temp2 = new File(filePath + img2);
		temp1.delete();
		temp2.delete();

		// Save as new image
		ImageIO.write(combined, "PNG", new File(path, "ReconTrace.png"));
	}

	private static void drawGuestTree(Tree tree, String fileName, boolean putLeafName, int[] size) throws IOException {

		// Getting Phylogeny object of the tree.
		final Phylogeny phy = getPhylogeny(tree, putLeafName, true);

		// Setting up a configuration object.
		final Configuration config = new Configuration();
		config.putDisplayColors(TreeColorSet.BACKGROUND, new Color(255, 255, 255));
		config.putDisplayColors(TreeColorSet.BRANCH, new Color(0, 0, 0));
		config.putDisplayColors(TreeColorSet.TAXONOMY, new Color(0, 0, 0));

		// Speciations: Red
		// Duplications: Green
		config.putDisplayColors(TreeColorSet.SPECIATION, new Color(255, 0, 0));
		config.putDisplayColors(TreeColorSet.DUPLICATION, new Color(0, 255, 0));
		config.putDisplayColors(TreeColorSet.DUPLICATION_OR_SPECATION, new Color(0, 0, 255));

		// Node Shape
		config.setDefaultNodeShape(NodeShape.RECTANGLE);
		config.setDefaultNodeFill(NodeFill.SOLID);
		config.setShowDefaultNodeShapesExternal(true);
		config.setShowDefaultNodeShapesInternal(true);
		config.setDefaultNodeShapeSize((short) 8);
		config.setBaseFontSize(14);

		// Branch related
		config.setColorizeBranches(false);
		config.setUseBranchesWidths(false);

		// Tree Style
		config.setPhylogenyGraphicsType(Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR);
		config.setDisplayNodeNames(true);
		config.setTaxonomyColorize(true);
		config.setDisplayTaxonomyCode(false);

		// Writing to a graphics file.
		AptxUtil.writePhylogenyToGraphicsFile(phy, new File(fileName), size[0], size[1], GraphicsExportType.PNG,
				config);
	}

	private static void drawHostTree(Tree tree, String fileName, boolean putLeafName, int[] size) throws IOException {

		// Getting Phylogeny object of the tree.
		final Phylogeny phy = getPhylogeny(tree, putLeafName, true);

		// Setting up a configuration object.
		final Configuration config = new Configuration();
		config.putDisplayColors(TreeColorSet.BACKGROUND, new Color(255, 255, 255));
		config.putDisplayColors(TreeColorSet.BRANCH, new Color(0, 0, 0));
		config.putDisplayColors(TreeColorSet.TAXONOMY, new Color(0, 0, 0));

		// Speciations = Red Duplications = Green Color
		config.putDisplayColors(TreeColorSet.SPECIATION, new Color(255, 0, 0));
		config.putDisplayColors(TreeColorSet.DUPLICATION, new Color(0, 255, 0));
		config.putDisplayColors(TreeColorSet.DUPLICATION_OR_SPECATION, new Color(0, 0, 255));

		// Node Shape
		config.setDefaultNodeShape(NodeShape.CIRCLE);
		config.setDefaultNodeFill(NodeFill.SOLID);
		config.setShowDefaultNodeShapesExternal(true);
		config.setShowDefaultNodeShapesInternal(true);
		config.setDefaultNodeShapeSize((short) 8);
		config.setBaseFontSize(14);

		// Branch related
		config.setColorizeBranches(false);
		config.setUseBranchesWidths(false);

		// Tree Style
		config.setPhylogenyGraphicsType(Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR);
		config.setDisplayNodeNames(true);
		config.setTaxonomyColorize(true);
		config.setDisplayTaxonomyCode(false);

		// Writing to a graphics file.
		AptxUtil.writePhylogenyToGraphicsFile(phy, new File(fileName), size[0], size[1], GraphicsExportType.PNG,
				config);

	}

	// Helper functions
	public static BufferedImage rotateImage(BufferedImage src, double degrees) {
		double radians = Math.toRadians(degrees);

		int srcWidth = src.getWidth();
		int srcHeight = src.getHeight();

		double sin = Math.abs(Math.sin(radians));
		double cos = Math.abs(Math.cos(radians));
		int newWidth = (int) Math.floor(srcWidth * cos + srcHeight * sin);
		int newHeight = (int) Math.floor(srcHeight * cos + srcWidth * sin);

		BufferedImage result = new BufferedImage(newWidth, newHeight, src.getType());
		Graphics2D g = result.createGraphics();
		g.translate((newWidth - srcWidth) / 2, (newHeight - srcHeight) / 2);
		g.rotate(radians, srcWidth / 2, srcHeight / 2);
		g.drawRenderedImage(src, null);

		return result;
	}

	public static void show(Tree tree, String fileName) throws IOException {

		boolean putLeafName = true;
		final Phylogeny phy = getPhylogeny(tree, putLeafName, false);

		// Displaying the newly created tree with Archaeopteryx.
		Archaeopteryx.createApplication(phy);
	}

	public static void writeTree(Tree tree, String fileName) throws IOException {

		int[] resolution = { 750, 750 };
		boolean putOnlyLeafName = false;
		boolean reconciledTree = false;

		// Getting Phylogeny object of the tree.
		final Phylogeny phy = getPhylogeny(tree, putOnlyLeafName, reconciledTree);

		// Setting up a configuration object.
		final Configuration config = new Configuration();
		config.putDisplayColors(TreeColorSet.BACKGROUND, new Color(255, 255, 255));
		config.putDisplayColors(TreeColorSet.BRANCH, new Color(0, 0, 0));
		config.putDisplayColors(TreeColorSet.TAXONOMY, new Color(0, 0, 0));

		// Speciations: Red
		config.putDisplayColors(TreeColorSet.SPECIATION, new Color(255, 0, 0));
		config.putDisplayColors(TreeColorSet.DUPLICATION, new Color(0, 255, 0));
		config.putDisplayColors(TreeColorSet.DUPLICATION_OR_SPECATION, new Color(0, 0, 255));

		// Node Shape
		config.setDefaultNodeShape(NodeShape.RECTANGLE);
		config.setDefaultNodeFill(NodeFill.SOLID);
		config.setShowDefaultNodeShapesExternal(true);
		config.setShowDefaultNodeShapesInternal(true);
		config.setDefaultNodeShapeSize((short) 6);
		config.setBaseFontSize(14);

		// Branch related
		config.setColorizeBranches(false);
		config.setUseBranchesWidths(false);

		// Tree Style
		config.setPhylogenyGraphicsType(Options.PHYLOGENY_GRAPHICS_TYPE.ROUNDED);
		config.setDisplayNodeNames(true);
		config.setDisplayInternalData(true);
		config.setTaxonomyColorize(true);
		config.setDisplayTaxonomyCode(false);

		// Writing to a graphics file.
		AptxUtil.writePhylogenyToGraphicsFile(phy, new File(fileName), resolution[0], resolution[1],
				GraphicsExportType.PNG, config);
	}

	// Creating a new rooted tree with given tree topology with its events.
	public static Phylogeny getPhylogeny(Tree tree, boolean putLeafName, boolean reconciled) throws IOException {

		final Phylogeny phy = new Phylogeny();

		final PhylogenyNode root;
		if (reconciled)
			root = constructReconciledTree(tree.root, null, putLeafName);
		else
			root = constructSimpleTree(tree.root, null, putLeafName);
		phy.setRoot(root);
		phy.setRooted(true);

		return phy;
	}

	private static PhylogenyNode constructSimpleTree(Node node, PhylogenyNode pnode, boolean putOnlyLeafName) {

		PhylogenyNode rootNode;
		if (!node.isLeaf()) {

			// First processing root node
			if (pnode == null) {
				rootNode = new PhylogenyNode();
				rootNode.setName(String.valueOf(node.id));
			} else {
				rootNode = pnode;
			}

			// Processing its children
			for (int i = 0; i < node.nchildren; i++) {
				Node child = node.children.get(i);
				PhylogenyNode d1 = new PhylogenyNode();
				d1.setName(String.valueOf(child.id));

				BranchData bd = new BranchData();

				bd.setBranchColor(new BranchColor(new Color(190, 190, 180)));
				bd.setBranchWidth(new BranchWidth(32.8));
				d1.setBranchData(bd);

				rootNode.addAsChild(constructSimpleTree(child, d1, putOnlyLeafName));
			}

		} else {
			final PhylogenyNode leaf = new PhylogenyNode();
			if (putOnlyLeafName)
				leaf.setName(node.name);
			else
				leaf.setName(String.valueOf(node.id) + " " + node.name);

			return leaf;
		}

		return rootNode;
	}

	/*
	 * private static PhylogenyNode constructSimpleTree(Node node,PhylogenyNode
	 * pnode, boolean putOnlyLeafName) {
	 * 
	 * PhylogenyNode rootNode; if (!node.isLeaf()) {
	 * 
	 * // First processing root node if(pnode==null) { rootNode = new
	 * PhylogenyNode(); rootNode.setName( String.valueOf(node.id) ); } else {
	 * rootNode = pnode; }
	 * 
	 * // Processing its children for(int i=0; i<node.nchildren ; i++) { Node child
	 * = node.children.get(i); PhylogenyNode d1 = new PhylogenyNode();
	 * d1.setName(String.valueOf(child.id));
	 * rootNode.addAsChild(constructSimpleTree(child,d1,putOnlyLeafName)); }
	 * 
	 * } else { final PhylogenyNode leaf = new PhylogenyNode(); if(putOnlyLeafName)
	 * leaf.setName(node.name); else leaf.setName(String.valueOf(node.id) + " " +
	 * node.name);
	 * 
	 * return leaf; }
	 * 
	 * return rootNode; }
	 */
	private static PhylogenyNode constructReconciledTree(Node node, PhylogenyNode pnode, boolean putLeafName) {

		PhylogenyNode rootNode;
		if (!node.isLeaf()) {

			// First processing root node
			if (pnode == null) {
				rootNode = new PhylogenyNode();
				rootNode.setName(String.valueOf(node.id));
			} else {
				rootNode = pnode;
			}

			// Root is a speciation or duplication
			// Event ev = getNodeEvent(node);
			// rootNode.getNodeData().setEvent(ev);

			// Processing its children
			for (int i = 0; i < node.nchildren; i++) {
				Node child = node.children.get(i);
				PhylogenyNode d1 = new PhylogenyNode();
				d1.setName(String.valueOf(child.id));
				rootNode.addAsChild(constructReconciledTree(child, d1, putLeafName));
			}

		} else {
			final PhylogenyNode leaf = new PhylogenyNode();
			if (putLeafName)
				leaf.setName(String.valueOf(node.id) + " " + node.name);
			else
				leaf.setName(String.valueOf(node.id));

			return leaf;
		}

		return rootNode;
	}

	/*
	 * private static Event getNodeEvent(Node node) { Event ev = new Event();
	 * 
	 * if( node.getEvent() == Reconciliation.EVENT_SPEC ) { ev.setSpeciations( 1 );
	 * ev.setDuplications( 0 ); } if( node.getEvent() == Reconciliation.EVENT_DUPL )
	 * { ev.setSpeciations( 0 ); ev.setDuplications( 1 ); } if( node.getEvent() ==
	 * Reconciliation.EVENT_LEAF ) { ev.setSpeciations( 0 ); ev.setDuplications( 0
	 * ); } return ev; }
	 */

	/*
	 * private static BranchData getBranchData(Node node) { BranchData bd = new
	 * BranchData(); bd.setBranchColor(new BranchColor(new Color( 170, 170, 170 )
	 * )); bd.setBranchWidth(new BranchWidth(65)); return bd; }
	 */
}
