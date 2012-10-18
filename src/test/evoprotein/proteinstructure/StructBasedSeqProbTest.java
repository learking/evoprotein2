package test.evoprotein.proteinstructure;

import test.beast.BEASTTestCase;
import evoprotein.evolution.datatype.MutableSequence;
import evoprotein.proteinstructure.InputStructure;
import evoprotein.proteinstructure.StructBasedSeqProb;
import evoprotein.proteinstructure.SolventAccessibility;
import evoprotein.proteinstructure.StructureEnv;
import junit.framework.TestCase;

public class StructBasedSeqProbTest extends TestCase {
	StructureEnv structEnv;
	SolventAccessibility solventAccessibility;
	InputStructure inputStructure;
	StructBasedSeqProb sequenceStructCompatibility;
	MutableSequence mutableSeq;
	
    @Override
    protected void setUp() throws Exception {
        super.setUp();
        solventAccessibility = new SolventAccessibility();
        solventAccessibility.initAndValidate();
        structEnv = new StructureEnv();
        structEnv.initAndValidate();
        inputStructure = new InputStructure();
        inputStructure.initByName("structureEnv", structEnv, "solventAccessibility", solventAccessibility);
        sequenceStructCompatibility = new StructBasedSeqProb();
        sequenceStructCompatibility.initByName("inputStructure", inputStructure);
    
    }
    
    public void testCalcSeqStructLogP() throws Exception{
        MutableSequence mutableSeq = new MutableSequence(9);
        // here codonUtil assumes 64 possible codons, while in reality, only 61 exist
        // should address this problem here and now
		int [] intSeq = new int [] {1,1,1,0,2,0,3,0,1};
		mutableSeq.setSequence(intSeq);
    	//double expectedLogP = 0;
    	System.out.println(sequenceStructCompatibility.calcSeqStructLogP(mutableSeq));
    	//assertEquals(expectedLogP, 1.0, BEASTTestCase.PRECISION);
    }
    
}
