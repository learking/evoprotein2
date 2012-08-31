package beast.evolution.likelihood;

import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.PathTree;

@Description("Calculates likelihood of PathTree using ... (to be filled in)")
public class PathTreeLikelihood extends Distribution {
	
    public Input<PathTree> m_pathTree = new Input<PathTree>("PathTree", "PathTree with sequence data in the all nodes", Validate.REQUIRED);
    public Input<SiteModel.Base> m_pSiteModel = new Input<SiteModel.Base>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);
	
    SubstitutionModel.Base m_substitutionModel;
    protected SiteModel.Base m_siteModel;
    
    @Override
    public void initAndValidate() throws Exception {
        m_siteModel = m_pSiteModel.get();
        m_substitutionModel = m_siteModel.m_pSubstModel.get();
        
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
