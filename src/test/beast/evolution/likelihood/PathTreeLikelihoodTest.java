/**
 * 
 */
package test.beast.evolution.likelihood;


import org.junit.Test;

import test.beast.evoprotein2TestCase;

import beast.evolution.alignment.Alignment;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.Tree;

/**
 * @author kuangyu
 *
 */
public class PathTreeLikelihoodTest extends evoprotein2TestCase {

	// weekend task 1:
	// create a path tree with internal states assigned and path mapped out
	// use this tree to test pathTreeLikelihood
	
	Alignment data;
	Tree tree;
	
    @Override
    protected void setUp() throws Exception {
        super.setUp();
        data = getDummyAlignment();
        tree = getTree(data);
    }
	
	@Test
	public void testLikelihoodCalculation() throws Exception {
		PathTree pathTree = new PathTree();
		pathTree.initByName("initial", tree, "alignment", data);
		/*
		int leftNr = pathTree.getRoot().getLeft().getNr();
		int rightNr = pathTree.getRoot().getRight().getNr();
		System.out.println("left nr:" + leftNr + " right Nr:" + rightNr + " Root height is:" + pathTree.getRoot().getHeight() + " 7 height:" + pathTree.getNode(7).getHeight());
		*/
		int [] dummySeq = {1,1,0,0,0};
		pathTree.setDummySeqInternalNodes(dummySeq);
		// pathTree.showSequences();
		
		
	}

}
