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
	
	public Input<PathTree> m_pathTree = new Input<PathTree>("pathtree", "pathtree on which this operation is performed", Validate.REQUIRED);
	
    @Override
    public void initAndValidate() throws Exception {
    }
	
	/*
	 * @see beast.core.Operator#proposal()
	 */
	@Override
	public double proposal() {
		// register this operator with input PathTree
		PathTree pathTree = m_pathTree.get(this);
		
		// Is there any way to utilize old path density without recalculation?
		// for example, if the proposal gets rejected
		double oldPathDensity, newPathDensity, fHastingsRatio;
		
		/*
		pathTree.setDummySeqInternalNodesAll();
		pathTree.showSequences();
	    System.out.println("--------------------------------------------------------------------------");
		*/
		
		// change the tree and set the tree to be dirty
		pathTree.setSomethingIsDirty(true);
		
		// MCMC to test our implementation:
		// case 1: reject
		
		// case 2: accept
		
		// calculate new path density
		
		// 1. Pupko: propose ancestral states
		
		// 2. Nielsen: propose substitution events
		
		// combine Pupko and Nielsen's partial density
		
		// calculate hastings ratio
		newPathDensity = 0.1;
		oldPathDensity = 0.2;
		fHastingsRatio = newPathDensity / oldPathDensity;
		// to make sure MCMC will reject everytime
		fHastingsRatio = Double.NEGATIVE_INFINITY;
		
		// to make sure MCMC will accept everytime
		// fHastingsRatio = Double.POSITIVE_INFINITY;
		
		return fHastingsRatio;
	}

}
