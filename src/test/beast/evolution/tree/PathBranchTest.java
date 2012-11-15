package test.beast.evolution.tree;


import java.util.ArrayList;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import evoprotein.evolution.datatype.MutableSequence;
import evoprotein.evolution.substitution.SubstitutionEvent;

import beast.evolution.tree.PathBranch;
import beast.evolution.tree.SeqPath;

import junit.framework.TestCase;

public class PathBranchTest extends TestCase {
	
	PathBranch pathBranch;
	List<SubstitutionEvent> list0;
	List<SubstitutionEvent> list1;
	List<SubstitutionEvent> list2;
	MutableSequence parentSeq, childSeq;
	
	@Before
	protected void setUp() throws Exception {
		super.setUp();
		pathBranch = new PathBranch(3,0,0);
		
		list0 = new ArrayList<SubstitutionEvent>();
		list0.add(new SubstitutionEvent(3,1,0.6));
		
		list1 = new ArrayList<SubstitutionEvent>();
		list1.add(new SubstitutionEvent(3,0,0.1));
		list1.add(new SubstitutionEvent(0,3,0.2));
		
		list2 = new ArrayList<SubstitutionEvent>();
		//list2.add(new SubstitutionEvent(2,3,0.1));
		//list2.add(new SubstitutionEvent(3,0,0.2));
		
		pathBranch.setMutationPath(0, list0);
		pathBranch.setMutationPath(1, list1);
		pathBranch.setMutationPath(2, list2);
		
		parentSeq = new MutableSequence(3);
		int [] intParentSeq = new int [] {3,3,0};
		parentSeq.setSequence(intParentSeq);
		
		childSeq = new MutableSequence(3);
		int [] intChildSeq = new int [] {1,3,0};
		childSeq.setSequence(intChildSeq);
	}

	@Test
	public void testGetSeqPath() throws Exception{
		SeqPath seqPath = pathBranch.getSeqPath(parentSeq, childSeq);
		System.out.println(seqPath.toString());
		System.err.println(seqPath.existStopCodon());
		System.out.println("=================");
		SeqPath codonSeqPath = pathBranch.getCodonSeqPath(1, parentSeq.getCodonSeq(0), childSeq.getCodonSeq(0));
		System.out.println(codonSeqPath.toString());
		System.err.println(codonSeqPath.existStopCodon());
	}

}
