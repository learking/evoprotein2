package evoprotein.evolution.datatype;

import java.util.Arrays;

import evoprotein.evolution.datatype.CodonUtil.Nucleotide;
import evoprotein.evolution.substitution.Substitution;

public class MutableSequence {
	
	static CodonUtil codonUtil = new CodonUtil();
	
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
	
	public int getNucleotide(int position){
		return intSequence[position];
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
	
	public MutableSequence getCodonSeq(int startSite) {
		MutableSequence codonSeq = new MutableSequence(3);
		//  +2 or +3 ? to be tested
		codonSeq.setSequence(new int [] {intSequence[startSite], intSequence[ startSite + 1 ], intSequence[ startSite + 2 ]});
		return codonSeq;
	}
	
	public void substitute(Substitution substitution) throws Exception{
		int site = substitution.getSite();
		int stateBeforeChange = substitution.getStateBeforeChange();
		if(stateBeforeChange == intSequence[site]){
			intSequence[site] = substitution.getStateAfterChange();
		}else{
			System.out.println("state before change:" + stateBeforeChange + "seq site: " + intSequence[site]);
			throw new Exception("state before substitution should be equal to unchanged seq at this site");
		}
	}
	
	public void mutate(int site, int newNucleotide) {
		intSequence[site] = newNucleotide;
	}
	
	// translate
	public int[] toCodonArray() {
		int codonArrayLength = intSequence.length / 3;
		int[] codonArray = new int[codonArrayLength];
		if (intSequence.length % 3 == 0) {
			for (int i = 0; i < intSequence.length - 2; i += 3) {
				codonArray[i/3] = codonUtil.translate(this, i);
			}
		}
		else {
			System.err.println("seqence is not multiple of 3!");
		}
		return codonArray;
	}

	public boolean existStopCodon(){
		boolean stopCodonFlag = false;
		int stopCodonPosition = -1;
		for (int startSite = 0; startSite < (intSequence.length - 2) ; startSite = startSite + 3) {

			int thisCodon = codonUtil.translate(this, startSite);
			if(thisCodon == -1){
				stopCodonFlag = true;
				stopCodonPosition = startSite;
				//System.err.println("Stop codon position:" + stopCodonPosition + "stop codon:" + thisCodon);
				break;
			}
		}
		return stopCodonFlag;
	}

	
	/*
	 * deep copy
	 */
	public MutableSequence copy(){
		MutableSequence mutableSequence = new MutableSequence(intSequence.length);
		mutableSequence.setSequence(intSequence);
		return mutableSequence;
	}
	
	@Override
	public boolean equals(Object aThat){
		if (!(aThat instanceof MutableSequence)) return false;
		
		MutableSequence that = (MutableSequence) aThat;
		
		if(Arrays.equals(that.getSequence(), intSequence)){
			return true;
		}else{
			return false;
		}
		
	}
	
	public String toString(){
		return Arrays.toString(intSequence);
	}
}
