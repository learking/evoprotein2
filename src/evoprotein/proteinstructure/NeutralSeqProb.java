/**
 * 
 */
package evoprotein.proteinstructure;

import evoprotein.evolution.datatype.MutableSequence;
import beast.core.Input;
import beast.core.Plugin;
import beast.core.Input.Validate;
import beast.evolution.substitutionmodel.Frequencies;

/**
 * @author kuangyu
 *
 */
public class NeutralSeqProb extends Plugin {
	// Note: not the frequencies that serve as input to the model used in the proposing stage
    public Input<Frequencies> frequenciesInput =
            new Input<Frequencies>("frequencies", "our substitution model equilibrium state frequencies", Validate.REQUIRED);
    
    double logP;
    
    Frequencies frequencies;
    double freqA;
    double freqC;
    double freqG;
    double freqT;
    double Y;
    
    @Override
    public void initAndValidate(){
    	frequencies = frequenciesInput.get();
    }
    
	public double calcNeutralSeqProb(MutableSequence seq) throws Exception{
		logP = 0;
		double[]  freqs = frequencies.getFreqs();
        freqA = freqs[0];
        freqC = freqs[1];
        freqG = freqs[2];
        freqT = freqs[3];
		Y = 1 - freqT * freqA * freqA - freqT * freqA * freqG - freqT * freqG * freqA ;

		int codonNumber = seq.getCodonNumber();
		int[] nucleoCounts = seq.getNucleoCounts();
		// a c g t
		logP = nucleoCounts[0]*Math.log(freqA) + nucleoCounts[1]*Math.log(freqC) + nucleoCounts[2]*Math.log(freqG) + nucleoCounts[3]*Math.log(freqT) - codonNumber*Math.log(Y);
		
		return getLogP();
	}
	
	public double getLogP(){
		return logP;
	}
	
}
