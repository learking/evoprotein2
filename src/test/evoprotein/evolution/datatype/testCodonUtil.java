package test.evoprotein.evolution.datatype;


import org.junit.Test;

import evoprotein.evolution.datatype.CodonUtil;
import evoprotein.evolution.datatype.MutableSequence;

import junit.framework.TestCase;

public class testCodonUtil extends TestCase {

	@Test
	public void testCodonTranslation() {
		MutableSequence mutableSeq = new MutableSequence(9);
		int [] intSeq = new int [] {1,1,1,0,2,0,3,3,3};
		mutableSeq.setSequence(intSeq);
		CodonUtil codonUtil = new CodonUtil();
		
		int firstCodonNr = codonUtil.translate(mutableSeq, 0);
		assertEquals(42, firstCodonNr);
		
		int secondCodonNr = codonUtil.translate(mutableSeq, 3);
		assertEquals(4, secondCodonNr);
		
		int thirdCodonNr = codonUtil.translate(mutableSeq, 6);
		assertEquals(63, thirdCodonNr);
	}

}
