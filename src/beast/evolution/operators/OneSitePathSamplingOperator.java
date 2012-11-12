/**
 * 
 */
package beast.evolution.operators;

import java.util.ArrayList;
import java.util.List;

import evoprotein.evolution.substitution.SubstitutionEvent;
import beast.core.Description;
import beast.evolution.tree.Node;
import beast.evolution.tree.PathBranch;
import beast.evolution.tree.PathTree;
import beast.util.Randomizer;

/**
 * @author kwang2
 *
 */
@Description("One site version of pathSamplingOperator, used in every MCMC round except initialization stage")
public class OneSitePathSamplingOperator extends PathSamplingOperator {

    int seqLength;
	
    // lambda parameter for nextExponential()
    final double lambda = 1.0;
    
    @Override
    public void accept(){
    	super.accept();
    	// when accept, update the value of oldPathLogDensity
    	oldPathLogDensity = newPathLogDensity;
    }
    
	@Override
	public double proposal() {
		// reset fHastingsRatio
		fHastingsRatio = 0;
		newPathLogDensity = 0;
		
		// register this operator with input PathTree
		PathTree pathTree = m_pathTree.get(this);
		
		if(oldPathLogDensity == Double.NEGATIVE_INFINITY){
			// first time: sample all sites across the tree
			//pathSampling(pathTree);
			oldPathLogDensity = newPathLogDensity;
			// let MCMC accept it
			fHastingsRatio = 999999999;
		}else{
			// after first time:
			// step 1: sample path for one site only
			
			// if encounter stop codon: reject directly
			
			// if no stop codon: calculate hastingRatio in normal fashion
			
			//pathSampling(pathTree);	
			// step 2: get hastingsRatio (assuming oldPathLogDensity is still valid (for now) !!!)
			fHastingsRatio = oldPathLogDensity - newPathLogDensity;
		}
		
		// change the tree and set the tree to be dirty
		pathTree.setSomethingIsDirty(true);
		return fHastingsRatio;
	}

	
	
	public void pathSampling(PathTree pathTree, int seqSite){
		int rootNr = pathTree.getRoot().getNr();
		int sudoRootNr = 0;
		for (Node childNode : pathTree.getRoot().getChildren()) {
			if(!childNode.isLeaf()) {
				sudoRootNr = childNode.getNr();
			}
		}
		
		// do it just for one site
		PupkoOneSite(pathTree, seqSite);
		
		// NielsenOneSite
		// don't need to traverse tree (we can work on m_branches directly, since now we have internal states already)
		for (int branchNr = 0; branchNr < pathTree.getBranches().size(); branchNr++) {
			if((branchNr != rootNr) && (branchNr != sudoRootNr)) {
				
				PathBranch thisBranch = pathTree.getBranch(branchNr);
				double thisBranchLength = pathTree.getNode(thisBranch.getEndNodeNr()).getLength();
				int parentNucleoState = pathTree.getSequences().get(thisBranch.getBeginNodeNr()).getSequence()[seqSite];
				int childNucleoState = pathTree.getSequences().get(thisBranch.getEndNodeNr()).getSequence()[seqSite];
				
				NielsenSampleOneBranchOneSite(thisBranch, seqSite, thisBranchLength, parentNucleoState, childNucleoState);
			}
		}
		
	}	
		
}
