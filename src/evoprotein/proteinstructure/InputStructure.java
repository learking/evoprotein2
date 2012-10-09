package evoprotein.proteinstructure;

import beast.core.Input;
import beast.core.Plugin;
import beast.core.Input.Validate;

public class InputStructure extends Plugin {
	
	public Input<StructureEnv> structureEnv = new Input<StructureEnv>("structureEnv", "Structure Environment", Validate.REQUIRED);
	
	int AAseqLength;
	
	double [] positionSpecificProb; 
	
	// need a sparse matrix here to store pairwise info, however, since the # of entries is much smaller than ...
	int [][] interactionTerm2EnvMap;
	
	// initiate and validate
	public void initAndValidate(){
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
		
		
		
	}
	
	// second order
	public void parseInteractionTerm2EnvMap(){
		
		// sparse matrix might be a better choice, but for now, I will use matrix
		
		
	}
	
	// 
	
}
