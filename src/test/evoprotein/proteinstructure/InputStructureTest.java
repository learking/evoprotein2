/**
 * 
 */
package test.evoprotein.proteinstructure;


import java.io.IOException;

import org.junit.Test;

import test.beast.BEASTTestCase;

import evoprotein.proteinstructure.InputStructure;
import evoprotein.proteinstructure.SolventAccessibility;
import evoprotein.proteinstructure.StructureEnv;

import junit.framework.TestCase;

/**
 * @author kuangyu
 *
 */
public class InputStructureTest extends TestCase {

	StructureEnv structEnv;
	SolventAccessibility solventAccessibility;
	InputStructure inputStructure;
	
    @Override
    protected void setUp() throws Exception {
        super.setUp();
        solventAccessibility = new SolventAccessibility();
        solventAccessibility.initAndValidate();
        structEnv = new StructureEnv();
        structEnv.initAndValidate();
        //String interactionTermsFile = "/home/kuangyu/workspace/evoprotein2/src/evoprotein/proteinstructure/interactionTerms.dat";
        inputStructure = new InputStructure();
        //inputStructure.initByName("structureEnv", structEnv, "solventAccessibility", solventAccessibility,"interactionTerms", interactionTermsFile);
        
    }
	
	@Test
	public void test() {
		assertEquals(structEnv.getStructEnvNum(), 10, BEASTTestCase.PRECISION);
	}

	@Test
	public void testWriteTerms() throws IOException {
		inputStructure.writeFirstOrder();
		System.out.println();
		inputStructure.writeInteraction2EnvMap();
	}
	
	
	/*
	@Test
	public void testParse() throws IOException{
			inputStructure.parseInteractionTerm2EnvMap();
	}
	*/

}
