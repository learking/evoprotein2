package evoprotein.proteinstructure;

import beast.core.Input;
import beast.core.Plugin;
import beast.core.Input.Validate;
import beast.evolution.tree.PathTree;
import evoprotein.evolution.datatype.MutableSequence;

public class StructBasedSeqProb extends Plugin {	
	// required input: inputStructure
	public Input<InputStructure> m_inputStructure = new Input<InputStructure>("inputStructure", "input Structure pre-calculated", Validate.REQUIRED);

	InputStructure inputStructure;
	
	public void initAndValidate(){
		inputStructure = m_inputStructure.get();
	}
	
	public double calcSeqStructLogP (MutableSequence seq) throws Exception {	
		
		double logP = 0;
		int [] codonArray = seq.toCodonArray();
		
		logP += calcFirstOrderLogP(codonArray);
		logP += calcInteractionLogP(codonArray);
		
		return logP;
	}
	
	public double calcFirstOrderLogP(int[] codonArray){
		double firstOrderLogP = 0;
		for (int i = 0 ; i < codonArray.length; i++) {
			firstOrderLogP += Math.log(inputStructure.getFirstOrderProb(i, codonArray[i]));
		}
		return firstOrderLogP;
	}
	
	public double calcInteractionLogP(int[] codonArray){
		double interactionLogP = 0;
		for (int m = 0; m < codonArray.length - 1 ; m++) {
			for (int n = m+1 ; n < codonArray.length ; n++) {
				interactionLogP += Math.log(inputStructure.getInteractionProb(m, n, codonArray[m], codonArray[n]));
			}
		}
		return interactionLogP;
	}
	
}
