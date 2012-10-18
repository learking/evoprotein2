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
	
	public int[] getNucleoCounts(){
		int[] nucleoCounts = new int[4];
		// a c g t
		for (int i = 0; i < intSequence.length ; i++) {
	        switch (intSequence[i]) {
	        case 0:  nucleoCounts[0]++;
	        		 break;
            case 1:  nucleoCounts[1]++;
                     break;
            case 2:  nucleoCounts[2]++;
                     break;
            case 3:  nucleoCounts[3]++;
                     break;
	        }
		}
		return nucleoCounts;
	}
	
	public int getCodonNumber() throws Exception{
		if(intSequence.length % 3 == 0){
			return intSequence.length / 3 ;
		}else{
			throw new Exception("Remainder is not ZERO!");
		}
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
