package test.beast.evolution.tree;


import org.junit.Before;
import org.junit.Test;

import beast.evolution.tree.SeqPath;

import evoprotein.evolution.datatype.MutableSequence;

import junit.framework.TestCase;

public class SeqPathTest extends TestCase {

	MutableSequence parentSeq, childSeq;
	
	@Before
	protected void setUp() throws Exception {
		super.setUp();
		parentSeq = new MutableSequence(3);
		int [] intParentSeq = new int [] {0,3,2};
		parentSeq.setSequence(intParentSeq);
		
		childSeq = new MutableSequence(3);
		int [] intChildSeq = new int [] {1,3,0};
		parentSeq.setSequence(intChildSeq);
	}

	@Test
	public void test() {
		//SeqPath seqPath =  new SeqPath(parentSeq, childSeq, substitutions);
	}

}
