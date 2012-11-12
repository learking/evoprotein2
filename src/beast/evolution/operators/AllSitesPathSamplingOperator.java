/**
 * 
 */
package beast.evolution.operators;

import java.util.ArrayList;
import java.util.List;

import evoprotein.evolution.datatype.MutableSequence;
import evoprotein.evolution.substitution.SubstitutionEvent;

import beast.core.Description;
import beast.evolution.tree.Node;
import beast.evolution.tree.PathBranch;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.SeqPath;
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
		
		// do it until internal seqs does not contain stop codon
		PupkoAllSites(pathTree);
		
		if(existStopCodonInternalNodes(pathTree)){
			System.out.println("there shouldn't be any stop codon any more!");
		}
		
		// NielsenOneSite
		// don't need to traverse tree (we can work on m_branches directly, since now we have internal states already)
		for (int branchNr = 0; branchNr < pathTree.getBranches().size(); branchNr++) {
			if((branchNr != rootNr) && (branchNr != sudoRootNr)) {
				NielsenSampleOneBranch(pathTree, branchNr);
			}
		}
		
		System.out.println("=====================================================================");
		
	}
	
	void PupkoAllSites(PathTree pathTree){
		// go codon by codon, make sure no stop codon appears
		for (int startSite = 0; startSite < (seqLength - 2) ; startSite = startSite + 3) {
			PupkoOneSite(pathTree, startSite);
			PupkoOneSite(pathTree, startSite + 1);
			PupkoOneSite(pathTree, startSite + 2);
			while (existStopCodonThisSiteInternalNodes(pathTree, startSite)) {
				PupkoOneSite(pathTree, startSite);
				PupkoOneSite(pathTree, startSite + 1);
				PupkoOneSite(pathTree, startSite + 2);
			}
		}
	}
	
	
	public void NielsenSampleOneBranch(PathTree pathTree, int branchNr) {

		PathBranch thisBranch = pathTree.getBranch(branchNr);
		
		int endNodeNr = thisBranch.getEndNodeNr();
		int beginNodeNr = thisBranch.getBeginNodeNr();
		MutableSequence childSeq = pathTree.getSequences().get(endNodeNr);
		MutableSequence parentSeq = pathTree.getSequences().get(beginNodeNr);
		
		double thisBranchLength = pathTree.getNode(thisBranch.getEndNodeNr()).getLength();
		
		// do it once, first
		for(int seqSite = 0; seqSite < seqLength; seqSite++){
			int parentNucleoState = pathTree.getSequences().get(thisBranch.getBeginNodeNr()).getSequence()[seqSite];
			int childNucleoState = pathTree.getSequences().get(thisBranch.getEndNodeNr()).getSequence()[seqSite];
			NielsenSampleOneBranchOneSite(thisBranch, seqSite, thisBranchLength, parentNucleoState, childNucleoState);
		}
		
		for (int startSite = 0; startSite < (seqLength - 2) ; startSite = startSite + 3) {
			for(int i = 0; i < 3; i++) {
				int currentSite = startSite + i;
				int parentNucleoState = parentSeq.getSequence()[currentSite];
				int childNucleoState = childSeq.getSequence()[currentSite];
				NielsenSampleOneBranchOneSite(thisBranch, currentSite, thisBranchLength, parentNucleoState, childNucleoState);
			}
			while (existStopCodonThisSiteInternalNodes(pathTree, startSite)) {
				for(int i = 0; i < 3; i++) {
					int currentSite = startSite + i;
					int parentNucleoState = parentSeq.getSequence()[currentSite];
					int childNucleoState = childSeq.getSequence()[currentSite];
					NielsenSampleOneBranchOneSite(thisBranch, currentSite, thisBranchLength, parentNucleoState, childNucleoState);
				}
			}
			
		}
		
	}
	
	public void NielsenSampleOneBranchOneSite(PathBranch thisBranch, int seqSite, double thisBranchLength, int parentNucleoState, int childNucleoState){
		
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

	}
	
	// checkers
	boolean existStopCodonInternalNodes(PathTree pathTree){
		boolean stopCodonFlag = false;
		for (Integer internalNodeIndex : internalNodesNr) {
			MutableSequence currentSeq = pathTree.getSequences().get(internalNodeIndex.intValue());
			if(currentSeq.existStopCodon()){
				stopCodonFlag = true;
				break;
			}
		}
		return stopCodonFlag;
	}
	
	
	//for debugging only
	public double getPathLogDensity() {
		return newPathLogDensity;
	}
	
	public void setPathLogDenstiyToZero() {
		newPathLogDensity = 0;
	}
}
