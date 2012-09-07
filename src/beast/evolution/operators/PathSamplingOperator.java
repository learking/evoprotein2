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
	PathTree m_pathTree;
	
    @Override
    public void initAndValidate() throws Exception {
    	m_pathTree = m_pPathTree.get();
    }
	
	/*
	 * @see beast.core.Operator#proposal()
	 */
	@Override
	public double proposal() {
		// TODO Auto-generated method stub
		
		return Double.NEGATIVE_INFINITY;
	}

}
