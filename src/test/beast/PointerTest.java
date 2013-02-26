package test.beast;


import org.junit.Test;

import evoprotein.evolution.datatype.MutableSequence;

import junit.framework.TestCase;

public class PointerTest extends TestCase {

	@Test
	public void test() {
		MutableSequence mutableSeq_a = new MutableSequence(9);
		int [] intSeq = new int [] {1,1,1,0,2,0,3,3,3};
		mutableSeq_a.setSequence(intSeq);
		
		MutableSequence mutableSeq_b = new MutableSequence(mutableSeq_a.getSequence().length);
		mutableSeq_b.setSequence(mutableSeq_a.getSequence());
		mutableSeq_b.mutate(0, 3);
		System.out.println(mutableSeq_a.toString());
		System.out.println(mutableSeq_b.toString());
	}

}
