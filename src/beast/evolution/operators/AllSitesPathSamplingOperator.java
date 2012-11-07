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
 * @author kuangyu
 *
 */

@Description("Used only during initialization stage of MCMC run")
public class AllSitesPathSamplingOperator extends PathSamplingOperator {
	
    double oldPathLogDensity = Double.NEGATIVE_INFINITY;
    
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
		int rootNr = pathTree.getRoot().getNr();
		int sudoRootNr = 0;
		for (Node childNode : pathTree.getRoot().getChildren()) {
			if(!childNode.isLeaf()) {
				sudoRootNr = childNode.getNr();
			}
		}
		
		for (int seqSite = 0; seqSite < seqLength; seqSite ++) {
			PupkoOneSite(pathTree, seqSite);
		}
		
		// NielsenOneSite
		// don't need to traverse tree (we can work on m_branches directly, since now we have internal states already)
		for (int branchNr = 0; branchNr < pathTree.getBranches().size(); branchNr++) {
			if((branchNr != rootNr) && (branchNr != sudoRootNr)) {
				NielsenSampleOneBranch(pathTree, branchNr);
			}
		}
	}
	
	public double NielsenSampleOneBranch(PathTree pathTree, int branchNr) {

		PathBranch thisBranch = pathTree.getBranch(branchNr);
		double thisBranchLength = pathTree.getNode(thisBranch.getEndNodeNr()).getLength();
		
		
		for(int seqSite = 0; seqSite < seqLength; seqSite++){
			int parentNucleoState = pathTree.getSequences().get(thisBranch.getBeginNodeNr()).getSequence()[seqSite];
			int childNucleoState = pathTree.getSequences().get(thisBranch.getEndNodeNr()).getSequence()[seqSite];
			NielsenSampleOneBranchOneSite(thisBranch, seqSite, thisBranchLength, parentNucleoState, childNucleoState);
		}
		
		/*
		System.out.println("branch:" + branchNr);
		int parentNucleoState = pathTree.getSequences().get(thisBranch.getBeginNodeNr()).getSequence()[0];
		int childNucleoState = pathTree.getSequences().get(thisBranch.getEndNodeNr()).getSequence()[0];
		NielsenSampleOneBranchOneSite(thisBranch, 0, thisBranchLength, parentNucleoState, childNucleoState);
		*/
		return 0;
	}
	
	public double NielsenSampleOneBranchOneSite(PathBranch thisBranch, int seqSite, double thisBranchLength, int parentNucleoState, int childNucleoState){
		
		final int childState = childNucleoState;
		final int parentState = parentNucleoState;
		final double totalTime = thisBranchLength;

		List<SubstitutionEvent> substitutionEvents = new ArrayList<SubstitutionEvent>();
		
		// CDFs for sampling another different nucleotide
		final double [][] differentCDFs = getDifferentNucleoCDF();
		
		// parameters used in the calculation
		int lastNucleotide = -1;
		double currentTime = 0;
		int currentState = parentState;
		int beginNucleotide = -1;
		int endNucleotide = -1;
		double timeInterval = 0;
		
		if(parentState == childState){
			while(lastNucleotide != childState){
				// initialize stuff
				substitutionEvents.clear();			
				currentTime = 0;
				currentState = parentState;
				
				while(currentTime < totalTime){
					timeInterval = Randomizer.nextExponential(lambda);
					if(currentTime + timeInterval >= totalTime){
						currentTime = totalTime;
						lastNucleotide = currentState;
					}else{
						currentTime += timeInterval;
						beginNucleotide = currentState;
						currentState = Randomizer.randomChoice(differentCDFs[currentState]);
						endNucleotide = currentState;
						// create new substitutionevent and add it to events
						// System.out.println("begin:" + beginNucleotide + " end:" + endNucleotide);
						substitutionEvents.add(new SubstitutionEvent(beginNucleotide, endNucleotide, timeInterval));
					}
				}
			}

			// copy substitutionEvents to ...
			if(substitutionEvents.size() != 0){
				thisBranch.setMutationPath(seqSite, substitutionEvents);
			}else{
				thisBranch.getMutationPath(seqSite).clear();
			}
			
		}
		else{
			while (lastNucleotide != childState) {
				// initialize stuff
				substitutionEvents.clear();			
				currentTime = 0;
				currentState = parentState;
				
				boolean firstSample = true;
				while (currentTime < totalTime) {
					if (firstSample) {
						timeInterval = sampleFirstSubstitutionTime(totalTime);
						currentTime = currentTime + timeInterval;
						
						beginNucleotide = currentState;
						currentState = Randomizer.randomChoice(differentCDFs[currentState]);
						endNucleotide = currentState;
						substitutionEvents.add(new SubstitutionEvent(beginNucleotide, endNucleotide, timeInterval));
						
						firstSample = false;
					} else {
						timeInterval = Randomizer.nextExponential(lambda);
						if (currentTime + timeInterval >= totalTime) {
							currentTime = totalTime;
							lastNucleotide = currentState;
						} else {
							currentTime += timeInterval;
							
							beginNucleotide = currentState;
							currentState = Randomizer.randomChoice(differentCDFs[currentState]);
							endNucleotide = currentState;
							substitutionEvents.add(new SubstitutionEvent(beginNucleotide, endNucleotide, timeInterval));
						}
					}
				}
			}
			
			// copy substitutionEvents to ...
			if(substitutionEvents.size() != 0){
				thisBranch.setMutationPath(seqSite, substitutionEvents);
			}else{
				System.err.print("Error: expecting at least one substitution at this site on this branch!");
			}
		}

		
		/*
		if(childState == parentState) {
			if(substitutionEvents.size() > 0) {
				System.out.println("begin: " + parentState + " child: " + childState);
				for (SubstitutionEvent substEvent : substitutionEvents) {
					System.out.println(substEvent.toString());
				}
			}
		}
		*/
		addToPathLogDensity(calculateOneSiteLogP(childState, parentState, totalTime, substitutionEvents));
		return 0;
	}
	
	
	//for debugging only
	public double getPathLogDensity() {
		return newPathLogDensity;
	}
	
	public void setPathLogDenstiyToZero() {
		newPathLogDensity = 0;
	}
}