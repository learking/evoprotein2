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
	
	// translate
	public int[] toCodonArray() throws Exception {
		int codonArrayLength = intSequence.length / 3;
		if (intSequence.length % 3 == 0) {
			int[] codonArray = new int[codonArrayLength];
			// translate each codon
			CodonUtil codonUtil = new CodonUtil();
			for (int i = 0; i < intSequence.length - 2; i += 3) {
				codonArray[i/3] = codonUtil.translate(this, i);
			}
			return codonArray;
		}
		else {
			throw new Exception("Remainder is not ZERO!");
		}
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
