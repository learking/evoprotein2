package evoprotein.proteinstructure;

import beast.core.Input;
import beast.core.Plugin;
import beast.core.Input.Validate;
import beast.evolution.tree.PathTree;
import evoprotein.evolution.datatype.MutableSequence;

public class SeqStructCompatibility extends Plugin {	
	// required input: inputStructure
	public Input<InputStructure> m_inputStructure = new Input<InputStructure>("inputStructure", "input Structure pre-calculated", Validate.REQUIRED);

	public double calcSeqStructProb (MutableSequence seq) {	
		
		return 1.0;
	}
	
}
