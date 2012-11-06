/**
 * 
 */
package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.PathTree;

/**
 * @author kwang2
 *
 */
@Description("One site version of pathSamplingOperator, used in every MCMC round except initialization stage")
public class OneSitePathSamplingOperator extends Operator {

	public Input<PathTree> m_pathTree = new Input<PathTree>("pathtree", "pathtree on which this operation is performed", Validate.REQUIRED);
    public Input<SiteModel.Base> m_pSiteModel = new Input<SiteModel.Base>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);
	
    SubstitutionModel.Base m_substitutionModel;
    double fHastingsRatio, newPathLogDensity;
	
    // lambda parameter for nextExponential()
    final double lambda = 1.0;
    
    @Override
    public void initAndValidate() throws Exception {
    	m_substitutionModel = m_pSiteModel.get().m_pSubstModel.get();
    }
    
    /*
    @Override
    public void accept(){
    	super.accept();
    	// when accept, update the value of oldPathLogDensity
    	oldPathLogDensity = newPathLogDensity;
    	
    	// for debugging
    	//estimateParameters();
    }
    */
    
	@Override
	public double proposal() {
		// TODO Auto-generated method stub
		return 0;
	}

	
	
}
