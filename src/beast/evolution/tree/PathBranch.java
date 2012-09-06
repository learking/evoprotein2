/**
 * 
 */
package beast.evolution.tree;

import java.util.ArrayList;
import java.util.List;

import evoprotein.evolution.substitution.SubstitutionEvent;

/**
 * @author kwang2
 *
 */
public class PathBranch {
	
	List<List<SubstitutionEvent>> m_MutationPaths = new ArrayList<List<SubstitutionEvent>>();
	
	public PathBranch() {}
	
	public PathBranch(int sequenceLength) {
		for(int i=0; i<sequenceLength; i++){
			m_MutationPaths.add(new ArrayList<SubstitutionEvent>());
		}
	}
	
	// getters
	public List<SubstitutionEvent> getMutationPath (int seqSite){
		return m_MutationPaths.get(seqSite);
	}
}
