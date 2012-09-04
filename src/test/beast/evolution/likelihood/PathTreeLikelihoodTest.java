/**
 * 
 */
package test.beast.evolution.likelihood;


import org.junit.Test;

import evoprotein.evolution.substitution.SubstitutionEvent;

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
		int nodeNr = 5;
		int leftNr = pathTree.getNode(nodeNr).getLeft().getNr();
		int rightNr = pathTree.getNode(nodeNr).getRight().getNr();
		System.out.println("left nr:" + leftNr + " LH:"+ pathTree.getNode(nodeNr).getLeft().getHeight() + " right Nr:" + rightNr + " RH:"+ pathTree.getNode(nodeNr).getRight().getHeight() + " height is:" + pathTree.getNode(nodeNr).getHeight());
		*/
		
		int [] dummySeq = {1,1,0,0,0};
		pathTree.setDummySeqInternalNodes(dummySeq);
		pathTree.showSequences();
		pathTree.setDummyPathBranch(2, 1, new SubstitutionEvent(1,3,0.25));
		pathTree.setDummyPathBranch(3, 1, new SubstitutionEvent(1,3,0.2));
		pathTree.setDummyPathBranch(4, 1, new SubstitutionEvent(1,3,0.2));
		
		pathTree.showOneSitePath(1);
		
	}

}
