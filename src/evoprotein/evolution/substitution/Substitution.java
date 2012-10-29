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
