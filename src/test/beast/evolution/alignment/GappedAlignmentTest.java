/**
 * 
 */
package test.beast.evolution.alignment;


import java.util.List;
import java.util.Set;

import org.junit.Test;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.GappedAlignment;
import beast.evolution.alignment.Sequence;

import junit.framework.TestCase;

/**
 * @author kuangyu
 *
 */
public class GappedAlignmentTest extends TestCase {

    static public GappedAlignment getGappedAlignment() throws Exception {
        Sequence human = new Sequence("human", "---ATCCCT------GGTTTT");
        Sequence chimp = new Sequence("chimp", "ACGTTA---ACGCTACTA---");

        GappedAlignment data = new GappedAlignment();
        data.initByName("sequence", human, "sequence", chimp,
                "dataType", "nucleotide"
        );
        return data;
    }
	
	@Test
	public void test() throws Exception {
		GappedAlignment data = getGappedAlignment();
		List<List<Integer>> seqData = data.getCounts();
		// it works! now confirmed that gap value is 16
		for (List<Integer> tmpSeq : seqData) {
			System.out.println(tmpSeq);
		}
		Set<Integer> delPositions = data.getDeletionPositions();
		for (Integer delPosition : delPositions) {
			System.out.print(delPosition + " ");
		}
	}

}
