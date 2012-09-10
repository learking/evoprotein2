/**
 * 
 */
package beast.evolution.operators;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.evolution.tree.PathTree;

/**
 * @author kuangyu
 *
 */


public class PathSamplingOperator extends Operator {
	
	public Input<PathTree> m_pPathTree = new Input<PathTree>("pathtree", "pathtree on which this operation is performed", Validate.REQUIRED);
	
    @Override
    public void initAndValidate() throws Exception {
    }
	
	/*
	 * @see beast.core.Operator#proposal()
	 */
	@Override
	public double proposal() {
		PathTree m_pathTree = m_pPathTree.get();
		
		// anyway to utilize old path density without recalculation?
		// for example, if the proposal gets rejected
		double oldPathDensity, newPathDensity, fHastingsRatio;
		
		// calculate new path density
		
		// 1. Pupko: propose ancestral states
		
		// 2. Nielsen: propose substitution events
		
		// combine Pupko and Nielsen's partial density
		
		// calculate hastings ratio
		newPathDensity = 0.1;
		oldPathDensity = 0.2;
		fHastingsRatio = newPathDensity / oldPathDensity;
		
		// 
		return fHastingsRatio;
	}

}
