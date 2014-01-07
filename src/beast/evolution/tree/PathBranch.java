/**
 * 
 */
package beast.evolution.tree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.core.Description;
import evoprotein.evolution.datatype.MutableSequence;
import evoprotein.evolution.substitution.Substitution;
import evoprotein.evolution.substitution.SubstitutionEvent;

/**
 * @author Kuangyu Wang
 *
 */
@Description("")
public class PathBranch {
	int beginNodeNr;
	int endNodeNr;

	List<List<SubstitutionEvent>> m_MutationPaths = new ArrayList<List<SubstitutionEvent>>();
	
	public PathBranch() {}
	
	public PathBranch(int sequenceLength, int beginNodeNumber, int endNodeNumber) {

		beginNodeNr = beginNodeNumber;
		endNodeNr = endNodeNumber;
		
		for(int i = 0; i < sequenceLength; i++){
			m_MutationPaths.add(new ArrayList<SubstitutionEvent>());
		}
	}
	
	// getters
	public List<SubstitutionEvent> getMutationPath (int seqSite){
		return m_MutationPaths.get(seqSite);
	}
	
	public int getEndNodeNr(){
		return endNodeNr;
	}
	
	public int getBeginNodeNr(){
		return beginNodeNr;
	}

	public SeqPath getSeqPath(MutableSequence parentSeq, MutableSequence childSeq) throws Exception{
		
		List<Substitution> substitutions = new ArrayList<Substitution>();
		
		// combine all mutations for each site 
		for (int site = 0; site < m_MutationPaths.size() ; site++) {
			int currentNumOfSubst = m_MutationPaths.get(site).size();
			if(currentNumOfSubst != 0){
				// within this branch, set heights
				double cumulativeHeight = 0;
				for (int substIndex = 0; substIndex < currentNumOfSubst; substIndex++ ) {
					cumulativeHeight += m_MutationPaths.get(site).get(substIndex).getTimeInterval();
					substitutions.add(new Substitution(site, m_MutationPaths.get(site).get(substIndex), cumulativeHeight));
				}
			}
		}
		
		// System.out.println("before sort:" + substitutions.toString());
		// sort mutations based on height
		Collections.sort(substitutions);
		// System.out.println("after sort:" + substitutions.toString());
		
		// create a list of mutable sequences
		// Given parent seq and substitutions, construct seq path
		SeqPath seqPath =  new SeqPath(parentSeq, childSeq, substitutions);
		
		//List<MutableSequence> seqPath = new ArrayList<MutableSequence>();
		return seqPath;
	}
	
	// focus 
	public SeqPath getCodonSeqPath(int siteNr, MutableSequence parentSeq, MutableSequence childSeq) throws Exception{
		List<Substitution> substitutions = new ArrayList<Substitution>();
		
		int startSite = siteNr - ( siteNr % 3 );
		
		//System.out.println("start site:" + startSite);
		
		
		// combine all mutations for the three sites belonging to this triplet 
		for (int site = startSite; site < (startSite + 3) ; site++) {
			int currentNumOfSubst = m_MutationPaths.get(site).size();
			if(currentNumOfSubst != 0){
				// within this branch, set heights
				double cumulativeHeight = 0;
				for (int substIndex = 0; substIndex < currentNumOfSubst; substIndex++ ) {
					cumulativeHeight += m_MutationPaths.get(site).get(substIndex).getTimeInterval();
					// site - startSite: since we are creating codonSeq below
					substitutions.add(new Substitution((site - startSite), m_MutationPaths.get(site).get(substIndex), cumulativeHeight));
				}
			}
		}
	
		Collections.sort(substitutions);
		//System.out.println(parentSeq.toString() + childSeq.toString());

		SeqPath codonSeqPath = new SeqPath(parentSeq.getCodonSeq(startSite) , childSeq.getCodonSeq(startSite) , substitutions);
		return codonSeqPath; 
		
	}
	
	// getter for testing purpose
	public int getTotalNumSubstitutions(){
		int totalNumSubstitutions = 0;
		for (int i = 0; i < m_MutationPaths.size(); i++) {
			totalNumSubstitutions += m_MutationPaths.get(i).size();
		}
		return totalNumSubstitutions;
	}
	
	// setters
	public void setMutationPath(int mutationPathIndex, List<SubstitutionEvent> newMutationPath) {
		m_MutationPaths.get(mutationPathIndex).clear();
		for(SubstitutionEvent tmpSubstitutionEvent : newMutationPath) {
			m_MutationPaths.get(mutationPathIndex).add(tmpSubstitutionEvent.copy());
		}
	}
	
	// deep copy
	public PathBranch copy(){
		PathBranch pathbranch = new PathBranch(m_MutationPaths.size(), beginNodeNr, endNodeNr);
		// if any branch has any substitution event, add them into pathbranch
		int mutationPathIndex = -1;
		for(List<SubstitutionEvent> tmpMutationPath: m_MutationPaths) {
			mutationPathIndex++;
			pathbranch.setMutationPath(mutationPathIndex, tmpMutationPath);
		}
		return pathbranch;
	}

	public void adjustSubstitutionTimes(double scaleFactor) {
		
		for (int site = 0; site < m_MutationPaths.size() ; site++) {
			int currentNumOfSubst = m_MutationPaths.get(site).size();
			if(currentNumOfSubst != 0){
				// within this branch, adjust interval length for each substitution
				for (int substIndex = 0; substIndex < currentNumOfSubst; substIndex++ ) {
					//substitutions.add(new Substitution(site, m_MutationPaths.get(site).get(substIndex), cumulativeHeight));
					m_MutationPaths.get(site).get(substIndex).adjustTimeInterval(scaleFactor);
				}
			}
		}

	}
	
	String mutationPath2String() {
		String mutationPathStr = "";
		//for each none empty path, record its position and substitution events (could be more than one)
		for(int i = 0; i < m_MutationPaths.size(); i++){
			//do this only when there is at least one substitution at this site
			if(m_MutationPaths.get(i).size() > 0){
				//convert each substitution event to string: site Nr, previousNucleo, currNucleo, time interval;
				String tmpPathStr = "";
				for(int j = 0; j < m_MutationPaths.get(i).size(); j++){
					//need to fix the toString() method of SubstitutionEvent
					tmpPathStr += "(" + Integer.toString(i) + "," + m_MutationPaths.get(i).get(j).toString() + ")";
				}
				mutationPathStr += tmpPathStr;
			}
		}
		return mutationPathStr;
	}
	
	public String toString() {
		String branchStr = "";
		//do the following only when it is a branch that has substitution event(s)
		if(getTotalNumSubstitutions()!=0){
			//need to put in all that are needed in order to reconstruct a PathBranch
			//part I: begin and end node
			String twoEndsStr = "(" + Integer.toString(getBeginNodeNr()) + "," + Integer.toString(getEndNodeNr()) + ")";
			//part II: mutationPaths
			String mutationPathStr = mutationPath2String();
			//join two parts
			branchStr = twoEndsStr + mutationPathStr + "|";
		}
		return branchStr;
	}
	
}
