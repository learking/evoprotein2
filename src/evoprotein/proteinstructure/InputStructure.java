package evoprotein.proteinstructure;

import beast.core.Input;
import beast.core.Plugin;
import beast.core.Input.Validate;

public class InputStructure extends Plugin {
	
	public Input<StructureEnv> structureEnv = new Input<StructureEnv>("structureEnv", "Structure Environment", Validate.REQUIRED);
	
	double [] positionSoventAccessbility; 
	
	// initiate and validate
	public void initAndValidate(){
		parseInputStructure();
	}
	
	// file parser
	public void parseInputStructure(){
		
		// set p(i|c_sol)
		
		
		// set p(m, n | c_sol, c_distance)
		
	}
	
	// 
	
}
