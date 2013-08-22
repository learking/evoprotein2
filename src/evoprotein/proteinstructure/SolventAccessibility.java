package evoprotein.proteinstructure;

import beast.core.Plugin;
import beast.util.Randomizer;

public class SolventAccessibility extends Plugin {
	
	double [][] solventCategories;
	
	public void initAndValidate(){
		parseSolventAccessibility();
	}
	
	public void parseSolventAccessibility(){
		// for now, mock up a table
		solventCategories = mockUpSolventCategories(10);
	}
	
	// getters
	
	public double[][] getSolventCategories(){
		return solventCategories;
	}
	
	public int getNumOfCategories() {
		return solventCategories.length;
	}
	
	public double getLogProb(int category, int codonType){
		return solventCategories[category][codonType];
	}
	
	// needs to be tested
	public double [][] mockUpSolventCategories(int numOfCategories){
		int rowNumber = numOfCategories; 
		double [][] mockupmatrix = new double [numOfCategories][61];
		// log version of mockupmatrix
		double [][] mockupLogMatrix = new double [numOfCategories][61];
		
		for (int rowNr = 0; rowNr < rowNumber; rowNr++) {
				mockupmatrix[rowNr] = new double[]{0.00983798,0.01745548,0.00222048,0.01443315,
						0.00844604,0.01498576,0.00190632,0.01239105,
						0.01064012,0.01887870,
						0.00469486,0.00833007,0.00688776,
						0.01592816,0.02826125,0.00359507,0.02336796,
						0.01367453,0.02426265,0.00308642,0.02006170,
						0.01722686,0.03056552,0.00388819,0.02527326,
						0.00760121,0.01348678,0.00171563,0.01115161,
						0.01574077,0.02792876,0.00355278,0.02309304,
						0.01351366,0.02397721,0.00305010,0.01982568,
						0.01702419,0.03020593,0.00384245,0.02497593,
						0.00751178,0.01332811,0.00169545,0.01102042,
						0.02525082,0.04480239,0.00569924,0.03704508,
						0.02167816,0.03846344,0.00489288,0.03180369,
						0.02730964,0.04845534,0.00616393,0.04006555,
						0.01205015,0.02138052,0.00271978,0.01767859};
		}
		
		// transform matrix to LogMatrix
		for (int rowNr = 0; rowNr < rowNumber; rowNr++) {
			for (int colNr = 0; colNr < 61; colNr++) {
				mockupLogMatrix[rowNr][colNr] = Math.log(mockupmatrix[rowNr][colNr]);
			}
		}
		
		return mockupLogMatrix;
	}
	
}
