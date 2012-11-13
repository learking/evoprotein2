/**
 * 
 */
package beast.evolution.operators;

import beast.core.Description;
import beast.evolution.tree.PathTree;

/**
 * @author kuangyu
 *
 */

@Description("Used only during initialization stage of MCMC run")
public class AllSitesPathSamplingOperator extends PathSamplingOperator {
    
    @Override
    public void accept(){
    	super.accept();
    	// when accept, update the value of oldPathLogDensity
    	oldPathLogDensity = newPathLogDensity;
    }
    
	/*
	 * @see beast.core.Operator#proposal()
	 */
	@Override
	public double proposal() {
		// reset fHastingsRatio
		fHastingsRatio = 0;
		newPathLogDensity = 0;
		
		// register this operator with input PathTree
		PathTree pathTree = m_pathTree.get(this);
		
		if(oldPathLogDensity == Double.NEGATIVE_INFINITY){
			// first time: calculate newPathLogDensity
			pathSampling(pathTree);
			oldPathLogDensity = newPathLogDensity;
			// let MCMC accept it
			fHastingsRatio = 999999999;
		}else{
			// after first time:
			// step 1: sample path
			pathSampling(pathTree);	
			// step 2: get hastingsRatio (assuming oldPathLogDensity is still valid (for now) !!!)
			fHastingsRatio = oldPathLogDensity - newPathLogDensity;
		}
		
		// change the tree and set the tree to be dirty
		pathTree.setSomethingIsDirty(true);
		return fHastingsRatio;
	}
	
	public void pathSampling(PathTree pathTree){
		
		// do it until internal seqs does not contain stop codon
		PupkoAllSites(pathTree);
		
		if(existStopCodonInternalNodes(pathTree)){
			System.out.println("there shouldn't be any stop codon any more!");
		}
		
		System.err.println("Pupko part done!");
		
		// NielsenOneSite
		// don't need to traverse tree (we can work on m_branches directly, since now we have internal states already)
		
		for (int branchNr = 0; branchNr < pathTree.getBranches().size(); branchNr++) {
			if((branchNr != rootNr) && (branchNr != sudoRootNr)) {
				try {
					NielsenSampleOneBranch(pathTree, branchNr);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		System.err.println("Nielsen part done!");
		
	}
	
	//for debugging only
	public double getPathLogDensity() {
		return newPathLogDensity;
	}
	
	public void setPathLogDenstiyToZero() {
		newPathLogDensity = 0;
	}
}
