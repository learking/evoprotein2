/**
 * 
 */
package test.beast.evolution.alignment;


import java.util.List;

import org.junit.Test;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;

import junit.framework.TestCase;

/**
 * @author kuangyu
 *
 */
public class NucleotideAlignmentTest extends TestCase {

    static public Alignment getAlignment() throws Exception {
        Sequence human = new Sequence("human", "AA-AACCC-CGGGGT-TTT");
        Sequence chimp = new Sequence("chimp", "ACGT-ACGTACG-TACG-T");

        Alignment data = new Alignment();
        data.initByName("sequence", human, "sequence", chimp,
                "dataType", "nucleotide"
        );
        return data;
    }
	
	@Test
	public void test() throws Exception {
		Alignment data = getAlignment();
		List<List<Integer>> seqData = data.getCounts();
		// it works! now confirmed that gap value is 16
		for (List<Integer> tmpSeq : seqData) {
			System.out.println(tmpSeq);
		}
	}

}
