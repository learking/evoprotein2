/**
 * 
 */
package beast.evolution.operators;

import beast.core.Description;
import beast.core.Operator;

/**
 * @author kwang2
 *
 */
@Description("One site version of pathSamplingOperator, used in every MCMC round except initialization stage")
public class OneSitePathSamplingOperator extends Operator {

	/* (non-Javadoc)
	 * @see beast.core.Operator#proposal()
	 */
	@Override
	public double proposal() {
		// TODO Auto-generated method stub
		return 0;
	}

}
