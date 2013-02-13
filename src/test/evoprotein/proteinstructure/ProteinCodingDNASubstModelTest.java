/**
 * 
 */
package test.evoprotein.proteinstructure;


import org.junit.Before;
import org.junit.Test;

import test.beast.BEASTTestCase;

import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.ProteinCodingDNASubstModel;
import evoprotein.evolution.datatype.MutableSequence;
import evoprotein.proteinstructure.InputStructure;
import evoprotein.proteinstructure.NeutralSeqProb;
import evoprotein.proteinstructure.SolventAccessibility;
import evoprotein.proteinstructure.StructBasedSeqProb;
import evoprotein.proteinstructure.StructureEnv;

import junit.framework.TestCase;

/**
 * @author kuangyu
 *
 */
public class ProteinCodingDNASubstModelTest extends TestCase {

	MutableSequence seqI, seqJ;
	RealParameter realParameter;
	Frequencies frequencies;
	
	SolventAccessibility solventAccessibility;
	StructureEnv structEnv;
	InputStructure inputStructure;
	
	NeutralSeqProb neutralSeqProb;
	StructBasedSeqProb sequenceStructCompatibility;
	ProteinCodingDNASubstModel ourModel;
	
	protected void setUp() throws Exception {
		super.setUp();
		
		Double[] freqs = new Double[] {0.1, 0.2, 0.3, 0.4};
		realParameter = new RealParameter(freqs);
		
		frequencies = new Frequencies();
		frequencies.initByName("frequencies", realParameter);
		
		// neutral model
		neutralSeqProb = new NeutralSeqProb();
		neutralSeqProb.initByName("frequencies", frequencies);
       
		seqI = new MutableSequence(9);
		int [] intSeqI = new int [] {1,1,1,0,2,0,3,0,1};
		seqI.setSequence(intSeqI);
		
		seqJ = new MutableSequence(9);
		int [] intSeqJ = new int [] {1,1,1,0,2,0,3,1,1};
		seqJ.setSequence(intSeqJ);
		
        solventAccessibility = new SolventAccessibility();
        solventAccessibility.initAndValidate();
        structEnv = new StructureEnv();
        structEnv.initAndValidate();
        inputStructure = new InputStructure();
        inputStructure.initByName("structureEnv", structEnv, "solventAccessibility", solventAccessibility);
        
        // struct based model
        sequenceStructCompatibility = new StructBasedSeqProb();
        sequenceStructCompatibility.initByName("inputStructure", inputStructure);
        
		// our Model
		ourModel = new ProteinCodingDNASubstModel();
		ourModel.initByName("frequencies", frequencies, "inputStructure", inputStructure);
	}
	
	@Test
	public void testGetSubstitutionRate() throws Exception {
		double seqI_logP = neutralSeqProb.calcNeutralSeqProb(seqI);
		double seqJ_logP = neutralSeqProb.calcNeutralSeqProb(seqJ);
		double neutral_expectedResult = seqJ_logP - seqI_logP;
		
		//double neutral_ourResult = Math.log(ourModel.getSubstitutionRate(seqI, seqJ));
		//System.out.println("I:" + seqI_logP + " J:" + seqJ_logP + " J-I:" + expectedResult);
		//System.out.println("ourresult:" + ourResult);
		//assertEquals(ourResult, expectedResult, BEASTTestCase.PRECISION);
		
		double seqI_structProb = sequenceStructCompatibility.calcSeqStructLogP(seqI);
		double seqJ_structProb = sequenceStructCompatibility.calcSeqStructLogP(seqJ);
		double struct_expectedResult = seqJ_structProb - seqI_structProb;
		
		//System.out.println("I:" + seqI_structProb + " J:" + seqJ_structProb + "J-I:" + struct_expectedResult);
		//System.out.println("our result:" + struct_ourResult);
		
		double logTAU_expected = struct_expectedResult - neutral_expectedResult;
		
		double substitutionRate_expected = 0.2 * logTAU_expected / (1 - 1 / Math.exp(logTAU_expected));
		
		double substitutionRate_ourResult = ourModel.getSubstitutionRate(seqI, seqJ, seqI.toCodonArray());
		
		assertEquals(substitutionRate_expected, substitutionRate_ourResult, BEASTTestCase.PRECISION);
	}

}
