/**
 * 
 */
package test.beast.substitutionmodel;


import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.InstantHKY;

import junit.framework.TestCase;

/**
 * @author kuangyu
 *
 */
public class InstantHKYTest extends TestCase {

    public interface Instance {
        Double[] getPi();

        Double getKappa();

        double[] getExpectedResult();
    }
	
    // beta = -1
    // k = 2
    
    protected Instance test0 = new Instance() {
        public Double[] getPi() {
            return new Double[]{0.25, 0.25, 0.25, 0.25};
        }

        public Double getKappa() {
            return 2.0;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    -1.0, 0.25, 0.5, 0.25,
                    0.25, -1.0, 0.25, 0.5,
                    0.5, 0.25, -1.0, 0.25,
                    0.25, 0.50, 0.25, -1.0
            };
        }
    }; 
    
    protected Instance test1 = new Instance() {
        public Double[] getPi() {
            return new Double[]{0.25, 0.25, 0.25, 0.25};
        }

        public Double getKappa() {
            return 2.0;
        }
        
        public double[] getExpectedResult() {
            return new double[]{
            		   0.933260716269620,   0.016901545023513,   0.032936193683354,   0.016901545023513,
            		   0.016901545023513,   0.933260716269620,   0.016901545023513,   0.032936193683354,
            		   0.032936193683354,   0.016901545023513,   0.933260716269620,   0.016901545023513,
            		   0.016901545023513,   0.032936193683354,   0.016901545023513,   0.933260716269620
            };
        }
    };
    
    Instance[] all = {test0};

    
	@Test
	public void testInstantHKY() throws Exception {
        for (Instance test : all) {

        	RealParameter f = new RealParameter(test.getPi());

        	Frequencies freqs = new Frequencies();
        	freqs.initByName("frequencies", f, "estimate", false);
        
        	InstantHKY instantHKY = new InstantHKY();
        	instantHKY.initByName("kappa", test.getKappa().toString(), "frequencies", freqs);

        	double[] mat = new double[4 * 4];
        	instantHKY.getInstantRateMatrix(mat);
        	final double[] result = test.getExpectedResult();
        
        	for (int k = 0; k < mat.length; ++k) {
        		assertEquals(mat[k], result[k], 1e-10);
        		System.out.println(k + " : " + (mat[k] - result[k]));
        	}
        }
	}

	@Test
	public void testTransitionProb() throws Exception{
    	RealParameter f = new RealParameter(test1.getPi());

    	Frequencies freqs = new Frequencies();
    	freqs.initByName("frequencies", f, "estimate", false);
    
    	InstantHKY instantHKY = new InstantHKY();
    	instantHKY.initByName("kappa", test1.getKappa().toString(), "frequencies", freqs);

    	double[] mat = new double[4 * 4];
    	instantHKY.getTransitionProbabilities(null, 0.07, 0, 1, mat);
    	final double[] result = test1.getExpectedResult();
    
    	for (int k = 0; k < mat.length; ++k) {
    		assertEquals(mat[k], result[k], 1e-10);
    		System.out.println(k + " : " + (mat[k] - result[k]));
    	}
	}
	
}
