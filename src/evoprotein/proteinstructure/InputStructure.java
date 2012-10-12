package evoprotein.proteinstructure;

import beast.core.Input;
import beast.core.Plugin;
import beast.core.Input.Validate;
import beast.util.Randomizer;

public class InputStructure extends Plugin {
	
	public Input<StructureEnv> structureEnv = new Input<StructureEnv>("structureEnv", "Structure Environment", Validate.REQUIRED);
	public Input<SolventAccessibility> solventAccessibility = new Input<SolventAccessibility>("solventAccessibility", "solvent accessibility categories", Validate.REQUIRED);
	
	int AAseqLength;
	
	int [] firstOrderTerms; 
	
	// need a sparse matrix here to store pairwise info, however, since the # of entries is much smaller than ...
	int [][] interactionTerm2EnvMap;
	
	// initiate and validate
	public void initAndValidate(){
		// need to initiate AA seq length
		AAseqLength = 100;

		parseInputStructure();
	}
	
	// file parser
	public void parseInputStructure(){
		
		// set p(i|c_sol)
		parseFirstOrderTerms();
		
		// set p(m, n | c_sol, c_distance)
		parseInteractionTerm2EnvMap();
	}
	
	// first order
	public void parseFirstOrderTerms(){
		// for now, mock up a list that logs which solvent category this position belongs to
		// length of the list should be equal to AAseqLength
		firstOrderTerms = mockUpFirstOrderTerms();
	}
	
	// second order
	public void parseInteractionTerm2EnvMap(){
		// sparse matrix might be a better choice, but for now, I will use matrix
		// for now, mock up a list that logs which environment this pair of AA under consideration belongs to
		interactionTerm2EnvMap = mockUpInteractionTerm2EnvMap();
	}
	
	// getters
	
	public double getFirstOrderProb(int codonPosition, int codonType) {
		int firstOrderTermCategory = firstOrderTerms[codonPosition];
		return solventAccessibility.get().getProb(firstOrderTermCategory, codonType); 
	}
	
	public double getInteractionProb(int firstCodonPosition, int secondCodonPosition, int firstCodonType, int secondCodonType){
		int structEnvNumber = interactionTerm2EnvMap[firstCodonPosition][secondCodonPosition];
		return structureEnv.get().getProb(structEnvNumber, firstCodonType, secondCodonType);
	}
	
	// for test only
	public int[] mockUpFirstOrderTerms(){
		// for categories available, randomly assign category to each position
		int [] firstOrder = new int [AAseqLength];
		double [] pdf = new double [solventAccessibility.get().getNumOfCategories()];
		
		double tempSum = 0;
		for (int i = 0; i < pdf.length; i++) {
			tempSum += 1.0 / (double) pdf.length;
			pdf[i] = tempSum;
		}
		
		for (int i = 0; i < firstOrder.length; i++) {
			firstOrder[i] = Randomizer.randomChoice(pdf);
		}
		
		return firstOrder;
	}
	
	public int[][] mockUpInteractionTerm2EnvMap(){
		int[][] interactionTerm2EnvMap = new int[AAseqLength][AAseqLength];

		int structEnvNum = structureEnv.get().getStructEnvNum();
		
		double tempSum = 0;
		double [] pdf = new double [structEnvNum];
		for (int i = 0; i < pdf.length; i++) {
			tempSum += 1.0 / (double) pdf.length;
			pdf[i] = tempSum;
		}
		
		for (int rowNr = 0 ; rowNr < interactionTerm2EnvMap.length; rowNr++) {
			for (int colNr = 0; colNr < interactionTerm2EnvMap[0].length; colNr++) {
				interactionTerm2EnvMap[rowNr][colNr] = Randomizer.randomChoice(pdf);
			}
		}
		
		return interactionTerm2EnvMap;
	}
	
}
