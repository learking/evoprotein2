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
		interactionTerm2EnvMap = mockUpInteractionTerm2EvnMap();
	}
	
	// for test only
	public int[] mockUpFirstOrderTerms(){
		// for categories available, randomly assign category to each position
		int [] firstOrder = new int [AAseqLength];
		double [] pdf = new double [solventAccessibility.get().numOfCategories()];
		
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
	
	public int[][] mockUpInteractionTerm2EvnMap(){
		int[][] interactionTerm2EvnMap = new int[AAseqLength][AAseqLength];

		// assign about 3N of all the entries in the matrix to sum up to "1"
		int nonZeroEntryNumber = AAseqLength * 3;
		
		
		
		return interactionTerm2EvnMap;
	}
	
}
