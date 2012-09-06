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
	
}
