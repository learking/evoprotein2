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

	@Test
	public void testToCodonArray() throws Exception {
		MutableSequence mutableSeq = new MutableSequence(9);
		int [] intSeq = new int [] {1,1,1,0,2,0,3,3,3};
		mutableSeq.setSequence(intSeq);
		int[] codonArray = mutableSeq.toCodonArray();
		int[] expectedArray = new int[3];
		expectedArray[0] = 42;
		expectedArray[1] = 4;
		expectedArray[2] = 63;
		for (int i = 0; i < codonArray.length; i++) {
			assertEquals(codonArray[i], expectedArray[i]);
		}
	}

}
