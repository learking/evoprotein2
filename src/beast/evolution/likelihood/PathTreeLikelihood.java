package beast.evolution.likelihood;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

import evoprotein.evolution.datatype.MutableSequence;
import evoprotein.evolution.substitution.SubstitutionEvent;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.InstantHKY;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.PathBranch;
import beast.evolution.tree.PathTree;

@Description("Calculates likelihood of PathTree using ... (to be filled in)")
public class PathTreeLikelihood extends Distribution {
	
    public Input<PathTree> m_pathTree = new Input<PathTree>("PathTree", "PathTree with sequence data in the all nodes", Validate.REQUIRED);
    public Input<SiteModel.Base> m_pSiteModel = new Input<SiteModel.Base>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);
	
    SubstitutionModel.Base m_substitutionModel;
    protected SiteModel.Base m_siteModel;
    // added for now
    InstantHKY m_substModel;
    
    @Override
    public void initAndValidate() throws Exception {
        m_siteModel = m_pSiteModel.get();
        m_substitutionModel = m_siteModel.m_pSubstModel.get();
        m_substModel = (InstantHKY) m_substitutionModel;
    }
    
    @Override
    public double calculateLogP() throws Exception{
    	System.out.println("ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc");
    	double oneSiteP = calculateOneSiteLogP(0);
    	logP = oneSiteP;
    	//return oneSiteP;
    	if(this.isDirtyCalculation()){
    		return logP;
    	}else{
    		return 100.0;
    	}
    }
    
    public double calculateOneSiteLogP(int seqSite) throws Exception{
		
    	double oneSiteLogP = 0.0;
		double[] instantMatrix = new double[4 * 4];
		m_substModel.getInstantRateMatrix(instantMatrix);
    	
		int sudoRootNr = 0;
		for (Node childNode : m_pathTree.get().getRoot().getChildren()) {
			if(!childNode.isLeaf()) {
				sudoRootNr = childNode.getNr();
			}
		}
		
    	for (int i=0 ; i< m_pathTree.get().getNodeCount(); i++){
    		if((!m_pathTree.get().getNode(i).isRoot()) && (m_pathTree.get().getNode(i).getNr()!=sudoRootNr)){
    			if(m_pathTree.get().getNode(i).getHeight() != m_pathTree.get().getRoot().getHeight()){
    				PathBranch currentBranch = m_pathTree.get().getBranch(i);

    				// sanity check
    				if (i != currentBranch.getEndNodeNr()){
    					throw new Exception("EndNode number does not match!");
    				}
    				
    				// get "end" and "begin" nucleotide
    				int endNodeNr = currentBranch.getEndNodeNr();
    				int beginNodeNr = currentBranch.getBeginNodeNr();
					MutableSequence endSeq = m_pathTree.get().getSequences().get(endNodeNr);
					MutableSequence beginSeq = m_pathTree.get().getSequences().get(beginNodeNr);
					int endNucleotide = endSeq.getSequence()[seqSite];
					int beginNucleotide = beginSeq.getSequence()[seqSite];
					int transitionOutEventCode;
					int transitionEventCode;
					double currentBranchLength = m_pathTree.get().getNode(endNodeNr).getLength();
					List<SubstitutionEvent> currentSubstitutionEvents = currentBranch.getMutationPath(seqSite);
    				
    				if(currentSubstitutionEvents.size() == 0){
    					// sanity check
    					// need to deal with first time calculation
    					
    					if(beginNucleotide != endNucleotide){
    						System.out.print("begin!=end");
    						//throw new Exception("begin and end nucleotide should be the same when there is no substitution along this branch!");
    					}
    					
    					// diagonal elements are already negative
    					transitionOutEventCode = beginNucleotide*4 + beginNucleotide;
    					oneSiteLogP += instantMatrix[transitionOutEventCode] * currentBranchLength;
    				}else{
    					// ...
    					double cumulativeHeight = 0.0;
    					int lastNucleotide = -1;
    					for (int substitutionIndex=0; substitutionIndex < currentSubstitutionEvents.size(); substitutionIndex ++){
    						
    						int previousNucleotide = currentSubstitutionEvents.get(substitutionIndex).getPreviousNucleotide();
    						int currentNucleotide = currentSubstitutionEvents.get(substitutionIndex).getCurrentNucleotide();
    						// determine the last nucleotide
    						if(substitutionIndex == (currentSubstitutionEvents.size()- 1)){
    							lastNucleotide = currentNucleotide;
    						}
    						
    						double currentTimeInterval = currentSubstitutionEvents.get(substitutionIndex).getTimeInterval();
    						cumulativeHeight += currentTimeInterval;
    						transitionOutEventCode = previousNucleotide * 4 + previousNucleotide;
    						transitionEventCode = previousNucleotide * 4 + currentNucleotide;
    						oneSiteLogP += Math.log(instantMatrix[transitionEventCode]) + instantMatrix[transitionOutEventCode] * currentTimeInterval;
    					}
    					
    					double unchangedHeight = currentBranchLength - cumulativeHeight;
    					// sanity check
    					
    					if(lastNucleotide != endNucleotide){
    						System.out.println("first:" + beginNucleotide + " end:" + endNucleotide + " last:" + lastNucleotide + " size:" + currentSubstitutionEvents.size());
    						System.out.println("last!=end");
    						throw new Exception("after last substituion, the nucleotide should be the same as the end of the branch! ");
    					}
    					
    					transitionOutEventCode = lastNucleotide*4 + lastNucleotide;
    					oneSiteLogP += instantMatrix[transitionOutEventCode] * unchangedHeight;
    				}
    			}
    		}
    	}
    	return oneSiteLogP;
    }
    
    /** CalculationNode methods **/
    
    @Override
    protected boolean requiresRecalculation() {
    	// if tree is dirty, recalculate
    	return m_pathTree.get().somethingIsDirty();
    }
    
    @Override
    public void store() {

        super.store();

    }

    @Override
    public void restore() {

        super.restore();

    }
    
	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub

	}

}
