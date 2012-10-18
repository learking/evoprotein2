/**
 * 
 */
package test.evoprotein.proteinstructure;

import test.beast.BEASTTestCase;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import evoprotein.evolution.datatype.MutableSequence;
import evoprotein.proteinstructure.NeutralSeqProb;
import junit.framework.TestCase;

/**
 * @author kuangyu
 *
 */
public class NeutralSeqProbTest extends  TestCase {

	MutableSequence mutableSeq;
	RealParameter realParameter;
	Frequencies frequencies;
	NeutralSeqProb neutralSeqProb;
	
	protected void setUp() throws Exception {
		super.setUp();
		
		Double[] freqs = new Double[] {0.1, 0.2, 0.3, 0.4};
		realParameter = new RealParameter(freqs);
		
		frequencies = new Frequencies();
		frequencies.initByName("frequencies", realParameter);
		
		neutralSeqProb = new NeutralSeqProb();
		neutralSeqProb.initByName("frequencies", frequencies);
       
		mutableSeq = new MutableSequence(9);
		int [] intSeq = new int [] {1,1,1,0,2,0,3,0,1};
		mutableSeq.setSequence(intSeq);
	}

	public void testCalcNeutralSeqProb() throws Exception {
		double [] freqs = frequencies.getFreqs();
		double freqA = freqs[0];
		double freqC = freqs[1];
		double freqG = freqs[2];
		double freqT = freqs[3];
		
		double expectedNeutralSeqProb = 3.0 *Math.log(freqA) + 4.0 *Math.log(freqC) + Math.log(freqG) + Math.log(freqT) - 3.0 * Math.log(1.0 - freqT * freqA * freqA - freqT * freqA * freqG - freqT * freqG * freqA);
		assertEquals(expectedNeutralSeqProb, neutralSeqProb.calcNeutralSeqProb(mutableSeq) , BEASTTestCase.PRECISION);
	}

}
