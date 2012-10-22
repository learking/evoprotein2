package beast.evolution.substitutionmodel;

import evoprotein.evolution.datatype.MutableSequence;
import evoprotein.proteinstructure.InputStructure;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.tree.Node;

public class ProteinCodingDNASubstModel extends CalculationNode {

	// determine inputs
	
	//
    //public Input<RealParameter> kappa = new Input<RealParameter>("kappa", "kappa parameter in our model", Validate.REQUIRED);
	
	// for use in neutral model
	public Input<Frequencies> frequenciesInput =
            new Input<Frequencies>("frequencies", "our substitution model equilibrium state frequencies", Validate.REQUIRED);
	// for use in struct based model
	public Input<InputStructure> m_inputStructure = new Input<InputStructure>("inputStructure", "input Structure pre-calculated", Validate.REQUIRED);

	// define variables here
    Frequencies frequencies;
    double freqA;
    double freqC;
    double freqG;
    double freqT;
    double Y;
	// initAndValidate
    @Override
    public void initAndValidate(){
    	frequencies = frequenciesInput.get();
    }
    
	public double getSubstitutionRate(MutableSequence seqI, MutableSequence seqJ) throws Exception{
		double substitutionRate = 0;

		// find where these two sequences differ (both location and value)
		int differPosition = getDifferPosition(seqI, seqJ);
		double tau = getStructBaseAndNeutralRatio(seqI, seqJ, differPosition);
		
		return substitutionRate;
	}
    
    int getDifferPosition(MutableSequence seqI, MutableSequence seqJ) throws Exception {
		int[] seq_i = seqI.getSequence();
		int[] seq_j = seqJ.getSequence();
    	int differPosition = -1;
    	int numOfDifferences = 0;
    	for (int nucleoPosition = 0; nucleoPosition < seq_i.length; nucleoPosition++) {
    		if(seq_i[nucleoPosition] != seq_j[nucleoPosition]){
    			numOfDifferences++;
    			differPosition = nucleoPosition;
    		}
    	}
    	if(numOfDifferences == 1){
    		return differPosition;
    	}
    	else{
    		throw new Exception("");
    	}
    }
    
    double getStructBaseAndNeutralRatio(MutableSequence seqI, MutableSequence seqJ, int differPosition){
    	double tau;
    	double neutralSeqProbRatio = frequencies.getFreqs()[seqJ.getSequence()[differPosition]] / frequencies.getFreqs()[seqJ.getSequence()[differPosition]];
    	double structBasedSeqProbRatio;
    	tau = neutralSeqProbRatio;
    	return tau;
    }
    
}
