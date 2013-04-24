package test.beast.evolution.operators;


import org.junit.Before;
import org.junit.Test;

import evoprotein.evolution.substitution.SubstitutionEvent;

import beast.core.State;
import beast.evolution.alignment.Alignment;
import beast.evolution.operators.BranchLengthOperator;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.Tree;

import test.beast.evoprotein2TestCase;

import junit.framework.TestCase;

public class BranchLengthOperatorTest extends evoprotein2TestCase{
	
	@Test
	public void testProposal() throws Exception {
		
		//Alignment data = getAlignmentWithNoStopCodon();
        Alignment data = getDummyAlignment();
        Tree tree = getTree(data);
        PathTree pathTree = new PathTree();
		pathTree.initByName("initial", tree, "alignment", data);
		
		// set seq paths
		pathTree.setDummyPathBranch(0, 1, new SubstitutionEvent(1,3,0.25));
		pathTree.setDummyPathBranch(1, 2, new SubstitutionEvent(1,3,0.25));
		pathTree.setDummyPathBranch(1, 3, new SubstitutionEvent(1,3,0.25));
		
		pathTree.setDummyPathBranch(3, 1, new SubstitutionEvent(1,3,0.25));
		pathTree.setDummyPathBranch(4, 1, new SubstitutionEvent(1,3,0.25));
		pathTree.setDummyPathBranch(4, 2, new SubstitutionEvent(1,3,0.25));
		pathTree.setDummyPathBranch(4, 3, new SubstitutionEvent(1,3,0.25));
		
		pathTree.setDummyPathBranch(5, 1, new SubstitutionEvent(1,3,0.25));		
		pathTree.setDummyPathBranch(6, 4, new SubstitutionEvent(1,3,0.25));
		pathTree.setDummyPathBranch(6, 5, new SubstitutionEvent(1,3,0.25));

		
		State state = new State();
		state.initByName("stateNode", pathTree);
		state.initialise();
		
		BranchLengthOperator branchOperator = new BranchLengthOperator();
		branchOperator.initByName("tree", pathTree);
		
		//for(int i=0; i < 100; i++){
			branchOperator.proposal();
			//fail("Not yet implemented");
		//}
	}

}
