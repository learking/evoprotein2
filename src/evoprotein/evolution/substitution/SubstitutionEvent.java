package evoprotein.evolution.substitution;

public class SubstitutionEvent {
	
	int previousNucleotide;
	int currentNucleotide;
	double timeInterval;
	
	public SubstitutionEvent(int beginNucleotide, int endNucleotide, double tInterval){
		previousNucleotide = beginNucleotide;
		currentNucleotide = endNucleotide;
		timeInterval = tInterval;
	}
	
	public String toString(){
		String outputString;
		outputString = Integer.toString(previousNucleotide) + "->" +Integer.toString(currentNucleotide) + " timeInterval: " + Double.toString(timeInterval);
		return outputString;
	}
	
	// Currently, we assume that if an event is not transversion, it must be transition
	public boolean isTransition(){
		if((previousNucleotide == 0 && currentNucleotide == 2)
				|| (previousNucleotide == 2 && currentNucleotide == 0)
				|| (previousNucleotide == 1 && currentNucleotide == 3)
				|| (previousNucleotide == 3 && currentNucleotide == 1)){
			return true;
		}else{
			return false;
		}
	}
	
	// getters
	public double getTimeInterval(){
		return timeInterval;
	}
	
	public int getCurrentNucleotide(){
		return currentNucleotide;
	}
	
	public int getPreviousNucleotide(){
		return previousNucleotide;
	}
	
	// copy
	public SubstitutionEvent copy() {
		SubstitutionEvent substitutionEvent = new SubstitutionEvent(previousNucleotide, currentNucleotide, timeInterval);
		return substitutionEvent;
	}
}
