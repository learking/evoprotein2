package evoprotein.evolution.datatype;

public class MutableSequence {
	
	int [] intSequence;
	
	public MutableSequence(int sequenceLength){
		intSequence = new int[sequenceLength];
	}
	
	// setters
	public void setSequence(int[] newSequence) throws Exception {
		// sanity check
		if(intSequence.length == newSequence.length){
			for(int i=0; i<newSequence.length; i++){
				intSequence[i] = newSequence[i];
			}
		}else{
			throw new Exception("new Sequence's length is different from old one!");
		}
	}
	
	// getters
	public int[] getSequence() {
		return intSequence;
	}
	
}
