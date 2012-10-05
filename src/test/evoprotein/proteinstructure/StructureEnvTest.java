package test.evoprotein.proteinstructure;


import org.junit.Test;

import test.beast.BEASTTestCase;

import evoprotein.proteinstructure.StructureEnv;

import junit.framework.TestCase;

public class StructureEnvTest extends TestCase {

	@Test
	public void testmockupMatrix() {
		StructureEnv structureEnv = new StructureEnv();
		double [][] mockupmatrix = structureEnv.mockUpMatrix();
		double matrixSum = 0;
		for (int i = 0; i < mockupmatrix.length; i++) {
			for (int j = 0; j < mockupmatrix[0].length; j++) {
				matrixSum += mockupmatrix[i][j];
			}
		}
		assertEquals(matrixSum, 1.0, BEASTTestCase.PRECISION);
	}

}
