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
}
