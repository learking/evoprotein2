/**
 * 
 */
package beast.evolution.tree;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;

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
		
		for(int i=0; i<sequenceLength; i++){
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
