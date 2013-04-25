package test.beast.evolution.operators;


import org.junit.Before;
import org.junit.Test;

import evoprotein.evolution.substitution.SubstitutionEvent;

import beast.core.State;
import beast.evolution.alignment.Alignment;
import beast.evolution.operators.BranchLengthOperator;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.SeqPath;
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
		
		// set seq for node 5, 6 and their parent
		pathTree.setDummySeq(5, new int[]{1,1,0,0,0,2});
		pathTree.setDummySeq(6, new int[]{1,3,0,0,1,1});
		pathTree.setDummySeq(7, new int[]{1,3,0,0,0,1});
		
		// set seq paths
		pathTree.setDummyPathBranch(0, 3, new SubstitutionEvent(2,0,0.01));
		pathTree.setDummyPathBranch(0, 3, new SubstitutionEvent(0,2,0.02));
		pathTree.setDummyPathBranch(0, 5, new SubstitutionEvent(2,1,0.03));
		pathTree.setDummyPathBranch(1, 5, new SubstitutionEvent(2,1,0.03));
		
		pathTree.setDummyPathBranch(5, 1, new SubstitutionEvent(3,1,0.1));
		pathTree.setDummyPathBranch(5, 5, new SubstitutionEvent(1,2,0.2));
			
		pathTree.setDummyPathBranch(4, 4, new SubstitutionEvent(1,0,0.2));
		pathTree.setDummyPathBranch(6, 4, new SubstitutionEvent(0,1,0.05));

		
		State state = new State();
		state.initByName("stateNode", pathTree);
		state.initialise();
		
		BranchLengthOperator branchOperator = new BranchLengthOperator();
		branchOperator.initByName("tree", pathTree);

		SeqPath path5 = pathTree.getBranch(5).getSeqPath(pathTree.getSequences().get(7), pathTree.getSequences().get(5));
		SeqPath path6 = pathTree.getBranch(6).getSeqPath(pathTree.getSequences().get(7), pathTree.getSequences().get(6));
		System.out.println(path5.toString());
		System.out.println(path6.toString());		
		//for(int i=0; i < 100; i++){
			branchOperator.proposal();
			//fail("Not yet implemented");
		//}
			SeqPath path5_new = pathTree.getBranch(5).getSeqPath(pathTree.getSequences().get(7), pathTree.getSequences().get(5));
			SeqPath path6_new = pathTree.getBranch(6).getSeqPath(pathTree.getSequences().get(7), pathTree.getSequences().get(6));
			System.out.println(path5_new.toString());
			System.out.println(path6_new.toString());	
	}

}
