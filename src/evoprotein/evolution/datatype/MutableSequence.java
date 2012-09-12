package evoprotein.evolution.datatype;

public class MutableSequence {
	
	int [] intSequence;
	
	public MutableSequence(int sequenceLength){
		intSequence = new int[sequenceLength];
	}
	
	// setters
	public void setSequence(int[] newSequence){
		// sanity check
		if(intSequence.length == newSequence.length){
			for(int i=0; i<newSequence.length; i++){
				intSequence[i] = newSequence[i];
			}
		}
	}
	
	// getters
	public int[] getSequence() {
		return intSequence;
	}
	
	/*
	 * deep copy
	 */
	public MutableSequence copy(){
		MutableSequence mutableSequence = new MutableSequence(intSequence.length);
		mutableSequence.setSequence(intSequence);
		return mutableSequence;
	}
}
