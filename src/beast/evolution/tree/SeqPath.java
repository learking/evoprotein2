package beast.evolution.tree;

import java.util.ArrayList;
import java.util.List;

import evoprotein.evolution.datatype.MutableSequence;
import evoprotein.evolution.substitution.Substitution;

public class SeqPath {

	List<MutableSequence> seqPath;
	
	double[] mutationTime;
	
	public SeqPath(MutableSequence parentSeq, MutableSequence childSeq , List<Substitution> substitutions) throws Exception{
		seqPath = new ArrayList<MutableSequence>();
		mutationTime = new double[substitutions.size()];
		
		MutableSequence tmpSeq = parentSeq.copy();
		for (int i = 0; i < substitutions.size(); i++) {
			mutationTime[i] = substitutions.get(i).getTime();
			tmpSeq.substitute(substitutions.get(i));
			seqPath.add(tmpSeq.copy());
		}
		
		if(!(seqPath.get(seqPath.size()-1).equals(childSeq))){
			throw new Exception("after last substitution, the seq should be the same as child seq");
		}
	}
	
	// getter
	public List<MutableSequence> getSeqs(){
		return seqPath;
	}
	
	public double[] getTimes(){
		return mutationTime;
	}
	
	// setter
	
	// toString
	public String toString(){
		String seqPathString = "";
		for (int i = 0; i < seqPath.size(); i++) {
			seqPathString += seqPath.get(i).toString() + " "+ mutationTime[i] + "\n";
		}
		return seqPathString;
	}
}
