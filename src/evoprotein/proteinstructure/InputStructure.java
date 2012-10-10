package evoprotein.proteinstructure;

import beast.core.Input;
import beast.core.Plugin;
import beast.core.Input.Validate;

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
	public int [] parseFirstOrderTerms(){
		// for now, mock up a list that logs which solvent category this position belongs to
		// length of the list should be equal to AAseqLength
		int [] firstOrder = new int [AAseqLength];
		double [] pdf = new double [solventAccessibility.get().numOfCategories()];
		
		
		return firstOrder;
	}
	
	// second order
	public void parseInteractionTerm2EnvMap(){
		// sparse matrix might be a better choice, but for now, I will use matrix
		// for now, mock up a list that logs which environment this pair of AA under consideration belongs to
		
		
	}
	
	// for test only
	public void mockUpFirstOrderTerms(){
		// for categories available, randomly assign category to each position
		
	}
	
	public void mockUpInteractionTerm2EvnMap(){
		
	}
	
}
