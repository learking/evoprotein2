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

@Description("This is a parent class, it contains utilities used by children pathSamplingOperators")
public abstract class PathSamplingOperator extends Operator {

	/* (non-Javadoc)
	 * @see beast.core.Operator#proposal()
	 */
	@Override
	public abstract double proposal();

}
