package test.evoprotein.proteinstructure;


import org.junit.Test;

import beast.util.Randomizer;

import test.beast.BEASTTestCase;

import evoprotein.proteinstructure.SolventAccessibility;

import junit.framework.TestCase;

public class SolventAccessibilityTest extends TestCase {

	@Test
	public void testMockUpSolventCategories() {
		SolventAccessibility solventAccessibility = new SolventAccessibility();
		double [][] mockupMatrix = solventAccessibility.mockUpSolventCategories(10);
		for (int rowNr = 0; rowNr < 10; rowNr++) {
			assertEquals(Randomizer.getTotal(mockupMatrix[rowNr]), 1.0, BEASTTestCase.PRECISION);
		}
	}

}
