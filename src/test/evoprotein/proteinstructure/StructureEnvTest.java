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

	@Test
	public void testMarginalProbCalculation(){
		StructureEnv structureEnv = new StructureEnv();
		double [][] mockupmatrix = structureEnv.mockUpMatrix();
		double [][] marginalMatrix = structureEnv.getMarginalProbMatrix(mockupmatrix);
		double rowTotal = 0;
		double colTotal = 0;
		for (int i = 0; i < marginalMatrix[0].length; i++) {
			rowTotal += marginalMatrix[0][i];
			colTotal += marginalMatrix[1][i];
		}
		assertEquals(colTotal, 1.0, BEASTTestCase.PRECISION);
	}
	
	@Test
	public void testLogProb(){
		StructureEnv structureEnv = new StructureEnv();
		structureEnv.parseStructureEnv();
		System.out.println(structureEnv.getLogProb(9, 43, 34));
	}
}
