/**
 * 
 */
package test.evoprotein.proteinstructure;


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
        inputStructure = new InputStructure();
        inputStructure.initByName("structureEnv", structEnv, "solventAccessibility", solventAccessibility);
    }
	
	@Test
	public void test() {
		assertEquals(structEnv.getStructEnvNum(), 10, BEASTTestCase.PRECISION);
	}

}
