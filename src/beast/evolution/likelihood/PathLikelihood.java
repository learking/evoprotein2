/**
 * 
 */
package beast.evolution.likelihood;

import java.util.List;
import java.util.Random;

import evoprotein.evolution.datatype.CodonUtil;
import evoprotein.evolution.datatype.MutableSequence;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.evolution.substitutionmodel.ProteinCodingDNASubstModel;
import beast.evolution.substitutionmodel.TwoStructSubstModel;
import beast.evolution.tree.PathBranch;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.SeqPath;
import beast.evolution.tree.Tree;

/**
 * @author kuangyu
 *
 */
@Description("Calculates likelihood of Path along the entire tree using (a rewrite of PathTreeLikelihood to fit our model)")
public class PathLikelihood extends Distribution {
	
	// inputs
    public Input<PathTree> m_pathTree = new Input<PathTree>("PathTree", "PathTree with sequence data in the all nodes", Validate.REQUIRED);
    
    // one struct
    //public Input<ProteinCodingDNASubstModel> m_ourModel = new Input<ProteinCodingDNASubstModel>("ourModel", "Our model that allows dependence among sites", Validate.REQUIRED);
    
    // two struct
    public Input<TwoStructSubstModel> m_ourModel = new Input<TwoStructSubstModel>("ourModel", "Our model that allows dependence among sites", Validate.REQUIRED);
    
	static CodonUtil codonUtil = new CodonUtil();
    
    PathTree pathTree;
    
    // one struct case
    //ProteinCodingDNASubstModel ourModel;
    
    // two structs case
    TwoStructSubstModel ourModel;
    
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
    	//System.out.println("Pathlikelihood:" + logP);
		return logP;
    }
	
    public double calcTotalPathLogP() throws Exception{
    	logP = 0;
    	// loop through all branches, except the one with length 0 (the beginning of the branch is rootSeq)
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
    	
    	// get root seq
    	int rootNr = pathTree.getRoot().getNr();
    	int [] rootCodonSeq = pathTree.getSequences().get(rootNr).toCodonArray();
    	// need to add rootSeq Stationary log prob back into the pathLikelihood
    	logP += ourModel.getRootSeqLogP(rootCodonSeq);
    	
    	return logP;
    }
    
    double calcPathLogP(PathBranch currentBranch, MutableSequence parentSeq, MutableSequence childSeq) throws Exception{
    	double pathLogP = 0;
    	
    	int childNodeNr = currentBranch.getEndNodeNr();
    	double currentBranchLength = pathTree.getNode(childNodeNr).getLength();
    	
    	SeqPath currentSeqPath = currentBranch.getSeqPath(parentSeq, childSeq);
    	/*
    	if(currentSeqPath.existStopCodon()){
    		throw new Exception("There cannot be any stop codon along the path at this point of time!");
    	}
    	*/
    	List<MutableSequence> currentSeqs = currentSeqPath.getSeqs();
    	double[] currentTimes = currentSeqPath.getTimes();
    	
    	if(currentTimes.length != 0){
    		// substitutions in between
    		// System.err.println("substitutions number:" + currentTimes.length);
    		for(int i = 0 ; i < currentTimes.length ; i++){
    			//System.out.println("start:" + System.currentTimeMillis());
    			if(i == 0){
    				// first substitution
    				int differPosition = getDifferPosition(parentSeq, currentSeqs.get(i));
    				int startSite = differPosition - (differPosition%3);
    				int differCodon = codonUtil.translate(currentSeqs.get(i), startSite);
    				
    				pathLogP += - ourModel.getSubstAwayRate(parentSeq) * currentTimes[i] + Math.log(ourModel.getSubstitutionRate(parentSeq, currentSeqs.get(i), parentSeq.toCodonArray(), differPosition, differCodon));
    			}else{
    				int differPosition = getDifferPosition(currentSeqs.get(i - 1), currentSeqs.get(i));
    				int startSite = differPosition - (differPosition%3);
    				int differCodon = codonUtil.translate(currentSeqs.get(i), startSite);
    				
    				pathLogP += - ourModel.getSubstAwayRate(currentSeqs.get(i - 1)) * (currentTimes[i] - currentTimes[i -1]) + Math.log(ourModel.getSubstitutionRate(currentSeqs.get(i - 1), currentSeqs.get(i), currentSeqs.get(i - 1).toCodonArray(), differPosition, differCodon));
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
    
    int getDifferPosition(MutableSequence seqI, MutableSequence seqJ) {
		int[] seq_i = seqI.getSequence();
		int[] seq_j = seqJ.getSequence();
    	int differPosition = -1;
    	//int numOfDifferences = 0;
    	for (int nucleoPosition = 0; nucleoPosition < seq_i.length; nucleoPosition++) {
    		if(seq_i[nucleoPosition] != seq_j[nucleoPosition]){
    			//numOfDifferences++;
    			differPosition = nucleoPosition;
    		}
    	}
    	return differPosition;
    	/*
    	if(numOfDifferences == 1){
    		return differPosition;
    	}
    	else{
    		throw new Exception("");
    	}
    	*/
    }
    
    @Override
    protected boolean requiresRecalculation() {
    	// check model's dirtiness first
    	
        if (ourModel.isDirtyCalculation()) {
            //m_nHasDirt = Tree.IS_DIRTY;
            return true;
        }
    	
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
