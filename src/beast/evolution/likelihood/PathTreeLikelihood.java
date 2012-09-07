package beast.evolution.likelihood;

import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.InstantHKY;
import beast.evolution.substitutionmodel.SubstitutionModel;
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
    	double oneSiteP = calculateOneSiteLogP(1);
    	return oneSiteP;
    }
    
    public double calculateOneSiteLogP(int seqSite) throws Exception{
    	double oneSiteLogP = 0.0;
		double[] instantMatrix = new double[4 * 4];
		m_substModel.getInstantRateMatrix(instantMatrix);
    	
    	for (int i=0 ; i< m_pathTree.get().getNodeCount(); i++){
    		if(!m_pathTree.get().getNode(i).isRoot()){
    			if(m_pathTree.get().getNode(i).getHeight() != m_pathTree.get().getRoot().getHeight()){
    				PathBranch currentBranch = m_pathTree.get().getBranch(i);
    				double currentHeight = 0.0;

    				// sanity check
    				if (i != currentBranch.getEndNodeNr()){
    					throw new Exception("EndNode number does not match!");
    				}
    				// TO-DO stuff
    				if(currentBranch.getMutationPath(seqSite).size() == 0){
    					// prob that nothing happened
    					
    					oneSiteLogP += 0.0;
    				}else{
    					// ...
    					
    					oneSiteLogP += 1.0;
    				}
    			}
    		}
    	}
    	return oneSiteLogP;
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
