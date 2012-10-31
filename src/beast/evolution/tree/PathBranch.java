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
		
		System.out.println("before sort:" + substitutions.toString());
		// sort mutations based on height
		Collections.sort(substitutions);
		System.out.println("after sort:" + substitutions.toString());
		
		// create a list of mutable sequences
		// Given parent seq and substitutions, construct seq path
		SeqPath seqPath =  new SeqPath(parentSeq, childSeq, substitutions);
		
		//List<MutableSequence> seqPath = new ArrayList<MutableSequence>();
		return seqPath;
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
	
}
