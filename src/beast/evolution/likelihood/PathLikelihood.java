/**
 * 
 */
package beast.evolution.likelihood;

import java.util.List;
import java.util.Random;

import evoprotein.evolution.datatype.MutableSequence;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.evolution.substitutionmodel.ProteinCodingDNASubstModel;
import beast.evolution.tree.PathBranch;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.SeqPath;

/**
 * @author kuangyu
 *
 */
@Description("Calculates likelihood of Path along the entire tree using (a rewrite of PathTreeLikelihood to fit our model)")
public class PathLikelihood extends Distribution {
	
	// inputs
    public Input<PathTree> m_pathTree = new Input<PathTree>("PathTree", "PathTree with sequence data in the all nodes", Validate.REQUIRED);
    public Input<ProteinCodingDNASubstModel> m_ourModel = new Input<ProteinCodingDNASubstModel>("ourModel", "Our model that allows dependence among sites", Validate.REQUIRED);
    
    PathTree pathTree;
    ProteinCodingDNASubstModel ourModel;
    
	// init
    @Override
    public void initAndValidate() throws Exception {
    	ourModel = m_ourModel.get();
    	pathTree = m_pathTree.get();
    }
    
    // a totally different calculation compared to PathTreeLikelihood
    @Override
    public double calculateLogP() throws Exception{
    	// if the first time calculating "oldLikelihood"
    	int rootNr = pathTree.getRoot().getNr();
    	int [] rootSeq = pathTree.getSequences().get(rootNr).getSequence();
    	boolean firstTimeCalculation;
    	int rootSeqTotal = 0;
    	for(int i = 0 ; i < rootSeq.length ; i++){
    		rootSeqTotal += rootSeq[0];
    	}
    	if(rootSeqTotal == 0){
    		firstTimeCalculation = true;
    	}
    	else{
    		firstTimeCalculation = false;
    	}
    	
    	if(firstTimeCalculation){
    		logP = -999999999;
    	}else{
        	logP = calcTotalPathLogP();
    	}
		return logP;
    }
	
    public double calcTotalPathLogP() throws Exception{
    	logP = 0;
    	// loop through all branches
    	for (int i=0 ; i< pathTree.getNodeCount(); i++){
			if(pathTree.getNode(i).getHeight() != pathTree.getRoot().getHeight()){
				PathBranch currentBranch = pathTree.getBranch(i);
				
				// sanity check
				if (i != currentBranch.getEndNodeNr()){
					throw new Exception("EndNode number does not match!");
				}
				
				int endNodeNr = currentBranch.getEndNodeNr();
				int beginNodeNr = currentBranch.getBeginNodeNr();
				MutableSequence childSeq = pathTree.getSequences().get(endNodeNr);
				MutableSequence parentSeq = pathTree.getSequences().get(beginNodeNr);
				
				// deal with each branch separately
				logP += calcPathLogP(currentBranch, parentSeq, childSeq);
				
			}
    	}
    	return logP;
    }
    
    double calcPathLogP(PathBranch currentBranch, MutableSequence parentSeq, MutableSequence childSeq) throws Exception{
    	double pathLogP = 0;
    	
    	int childNodeNr = currentBranch.getEndNodeNr();
    	double currentBranchLength = pathTree.getNode(childNodeNr).getLength();
    	
    	SeqPath currentSeqPath = currentBranch.getSeqPath(parentSeq, childSeq);
    	if(currentSeqPath.existStopCodon()){
    		throw new Exception("There cannot be any stop codon along the path at this point of time!");
    	}
    	List<MutableSequence> currentSeqs = currentSeqPath.getSeqs();
    	double[] currentTimes = currentSeqPath.getTimes();
    	
    	if(currentTimes.length != 0){
    		// substitutions in between
    		// System.err.println("substitutions number:" + currentTimes.length);
    		for(int i = 0 ; i < currentTimes.length ; i++){
    			//System.out.println("start:" + System.currentTimeMillis());
    			if(i == 0){
    				// first substitution
    				pathLogP += - ourModel.getSubstAwayRate(parentSeq) * currentTimes[i] + Math.log(ourModel.getSubstitutionRate(parentSeq, currentSeqs.get(i)));
    			}else{
    				pathLogP += - ourModel.getSubstAwayRate(currentSeqs.get(i - 1)) * (currentTimes[i] - currentTimes[i -1]) + Math.log(ourModel.getSubstitutionRate(currentSeqs.get(i - 1), currentSeqs.get(i)));
    			}
    			//System.out.println("end:" + System.currentTimeMillis());
    		}
    		
    		// last substitution
    		if(currentSeqs.get(currentSeqs.size()-1).equals(childSeq)){
    		pathLogP += - ourModel.getSubstAwayRate(currentSeqs.get(currentSeqs.size()-1)) * (currentBranchLength - currentTimes[currentTimes.length - 1]);
    		}else{
    			throw new Exception("during last interval, there should be no substitution");
    		}
    		
    	}else{
    		// no change in this branch
    		System.err.println("branch without substitution exists");
    		if(parentSeq.equals(childSeq)){
    			pathLogP += - ourModel.getSubstAwayRate(parentSeq) * currentBranchLength;
    		}else{
    			throw new Exception("begin and end seq differ, while no substituiton found");
    		}
    	}
		//System.err.println("one branch done!");
    	return pathLogP;
    }
    
    @Override
    protected boolean requiresRecalculation() {
    	// if tree is dirty, recalculate
    	return pathTree.somethingIsDirty();
    }
    
	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see beast.core.Distribution#getConditions()
	 */
	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see beast.core.Distribution#sample(beast.core.State, java.util.Random)
	 */
	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub

	}

}
