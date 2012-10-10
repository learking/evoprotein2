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
	
	public int numOfCategories() {
		return solventCategories.length;
	}
	
	// needs to be tested
	public double [][] mockUpSolventCategories(int numOfCategories){
		int rowNumber = numOfCategories; 
		double [][] mockupmatrix = new double [numOfCategories][61];
		
		for (int rowNr = 0; rowNr < rowNumber; rowNr++) {
			for (int colNr = 0; colNr < 61; colNr++) {
				mockupmatrix[rowNr][colNr] = Randomizer.nextDouble();
			}
		}
		
		// normalize each column
		for (int rowNr = 0; rowNr < rowNumber; rowNr++) {
			mockupmatrix[rowNr] = Randomizer.getNormalized(mockupmatrix[rowNr]);
		}
		
		return mockupmatrix;
	}
	
}
