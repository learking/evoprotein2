package evoprotein.proteinstructure;

import beast.core.Input;
import beast.core.Plugin;
import beast.core.Input.Validate;
import beast.evolution.tree.PathTree;
import evoprotein.evolution.datatype.MutableSequence;

public class SeqStructCompatibility extends Plugin {
	
	// required input: structureEnv
	public Input<StructureEnv> m_pathTree = new Input<StructureEnv>("structenv", "Structure Env pre-calculated", Validate.REQUIRED);
	// required input: inputStructure
	
	public double calcSeqStructProb (MutableSequence seq) {	
		
		
		return 1.0;
	}
	
}
