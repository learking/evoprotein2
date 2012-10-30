/**
 * 
 */
package test.evoprotein.evolution.datatype;


import org.junit.Test;

import evoprotein.evolution.datatype.MutableSequence;

import junit.framework.TestCase;

/**
 * @author kwang2
 *
 */
public class MutableSequenceTest extends TestCase {
	
	MutableSequence mutableSeq;
	
	@Override
	protected void setUp() {
		mutableSeq = new MutableSequence(9);
		int [] intSeq = new int [] {1,1,1,0,2,0,3,3,3};
		mutableSeq.setSequence(intSeq);
	}
	
	@Test
	public void testToCodonArray() throws Exception {
		int[] codonArray = mutableSeq.toCodonArray();
		int[] expectedArray = new int[3];
		expectedArray[0] = 42;
		expectedArray[1] = 4;
		expectedArray[2] = 60;
		for (int i = 0; i < codonArray.length; i++) {
			assertEquals(codonArray[i], expectedArray[i]);
		}
	}

	@Test
	public void testNucleoCounts(){
		int [] expectedCounts = new int [] {2, 3, 1, 3};
		int [] realCounts = mutableSeq.getNucleoCounts();
		for (int i = 0; i < expectedCounts.length; i++ ) {
			assertEquals(expectedCounts[i], realCounts[i]);
		}
	}
	
	@Test
	public void testEquals(){
		MutableSequence anotherSeq = new MutableSequence(9);
		int [] intAnotherSeq = new int [] {1,1,1,0,2,0,3,3,3};
		anotherSeq.setSequence(intAnotherSeq);
		assertEquals(true, anotherSeq.equals(mutableSeq));
	}
	
}
