package evoprotein.evolution.substitution;

public class Substitution  implements Comparable<Substitution> {
	
	int site;
	int stateBeforeChange;
	int stateAfterChange;
	double time;
	
	public Substitution(int substSite, SubstitutionEvent substitutionEvent, double substitutionTime){
		site = substSite;
		stateBeforeChange = substitutionEvent.getPreviousNucleotide();
		stateAfterChange = substitutionEvent.getCurrentNucleotide();
		time = substitutionTime;
	}

	// getter
	public double getTime(){
		return time;
	}
	
	public int getSite(){
		return site;
	}
	
	public int getStateBeforeChange(){
		return stateBeforeChange;
	}
	
	public int getStateAfterChange(){
		return stateAfterChange;
	}	
	
	// toString
	public String toString(){
		String outputString;
		outputString = Integer.toString(site) + "||" + Integer.toString(stateBeforeChange) + "->" +Integer.toString(stateAfterChange) + " time: " + Double.toString(time);
		return outputString;
	}
	
	// Comparable interface implementation
	@Override
	public int compareTo(Substitution n) {
		if(time < n.getTime()){
			return -1;
		}else if (time > n.getTime()) {
			return  1;
		}else{
			return 0;
		}
	}
	
	
}
