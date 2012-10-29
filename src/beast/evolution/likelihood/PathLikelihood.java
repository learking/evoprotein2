/**
 * 
 */
package beast.evolution.likelihood;

import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.evolution.substitutionmodel.ProteinCodingDNASubstModel;
import beast.evolution.tree.PathTree;

/**
 * @author kuangyu
 *
 */
@Description("Calculates likelihood of Path along the entire tree using (a rewrite of PathTreeLikelihood to fit our model)")
public class PathLikelihood extends Distribution {
	
	// inputs
    public Input<PathTree> m_pathTree = new Input<PathTree>("PathTree", "PathTree with sequence data in the all nodes", Validate.REQUIRED);
    public Input<ProteinCodingDNASubstModel> m_ourModel = new Input<ProteinCodingDNASubstModel>("ourModel", "Our model that allows dependence among sites", Validate.REQUIRED);
    
    ProteinCodingDNASubstModel ourModel;
    
	// init
    @Override
    public void initAndValidate() throws Exception {
    	ourModel = m_ourModel.get();
    }
    
    // a totally different calculation compared to PathTreeLikelihood
    @Override
    public double calculateLogP() throws Exception{
    	// if the first time calculating "oldLikelihood"
    	int rootNr = m_pathTree.get().getRoot().getNr();
    	int [] rootSeq = m_pathTree.get().getSequences().get(rootNr).getSequence();
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
        	logP = calcPathLogP();
    	}
		return logP;
    }
	
    public double calcPathLogP(){
    	
    	// combine substitutions among all sites into a single sequence path
    	
    	// loop through all branches
    	
    	// deal with each branch separately
    	
    	// add them all up
    	
    	return 1.0;
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
