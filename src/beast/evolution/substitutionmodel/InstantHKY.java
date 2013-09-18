/**
 * 
 */
package beast.evolution.substitutionmodel;

/**
 * @author kuangyu
 *
 */
public class InstantHKY extends HKY {

	public void getInstantRateMatrix(double [] matrix){
		
		// need a condition here
		setupMatrix();
		
		double k = kappaInput.get().getValue();
		double freqR = freqA + freqG;
		double freqY = freqC + freqT;
		
		// Row 1: A
		matrix[0] = k * beta * freqG + beta * freqY;
		matrix[1] = - beta * freqC;
		matrix[2] = - k * beta * freqG;
		matrix[3] = - beta * freqT;
		// Row 2: C
		matrix[4] = - beta * freqA;
		matrix[5] = k * beta * freqT + beta * freqR;
		matrix[6] = - beta * freqG;
		matrix[7] = - k* beta * freqT;
		// Row 3: G
		matrix[8] = - k * beta * freqA;
		matrix[9] = - beta * freqC;
		matrix[10] = k * beta * freqA + beta*freqY;
		matrix[11] = - beta * freqT;
		// Row 4: T
		matrix[12] = - beta * freqA;
		matrix[13] = - k * beta * freqC;
		matrix[14] = - beta * freqG;
		matrix[15] = k * beta * freqC + beta * freqR;
	}
	
}
