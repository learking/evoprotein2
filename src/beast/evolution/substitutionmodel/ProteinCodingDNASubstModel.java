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
    InputStructure inputStructure;
    double freqA;
    double freqC;
    double freqG;
    double freqT;
    double Y;
	// initAndValidate
    @Override
    public void initAndValidate(){
    	frequencies = frequenciesInput.get();
    	inputStructure = m_inputStructure.get();
    }
    
	public double getSubstitutionRate(MutableSequence seqI, MutableSequence seqJ) throws Exception{
		double substitutionRate = 0;

		// find where these two sequences differ (both location and value)
		int differPosition = getDifferPosition(seqI, seqJ);
		double logTAU = getLogTAU(seqI, seqJ, differPosition);
		substitutionRate = logTAU / (1 - 1/Math.exp(logTAU));
		if(isTransition(seqI.getNucleotide(differPosition), seqJ.getNucleotide(differPosition))){
			// should be
			//substitutionRate = substitutionRate * frequencies.getFreqs()[seqJ.getNucleotide(differPosition)];
			substitutionRate = substitutionRate * frequencies.getFreqs()[seqJ.getNucleotide(differPosition)];
		}else{
			substitutionRate = substitutionRate * frequencies.getFreqs()[seqJ.getNucleotide(differPosition)];
		}
		
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
    
    boolean isTransition(int firstNucleotide, int secondNucleotide){
    	int sum = firstNucleotide + secondNucleotide;
    	if(sum == 2 || sum == 4){
    		return true;
    	}else{
    		return false;
    	}
    }
    
    double getLogTAU(MutableSequence seqI, MutableSequence seqJ, int differPosition) throws Exception{
    	double logTAU;
    	double neutralSeqProbRatio = Math.log(frequencies.getFreqs()[seqJ.getSequence()[differPosition]] / frequencies.getFreqs()[seqI.getSequence()[differPosition]]);
    	double structBasedSeqProbRatio = getStructBasedSeqProbRatio(seqI, seqJ, differPosition);
    	logTAU = structBasedSeqProbRatio - neutralSeqProbRatio;
    	return logTAU;
    }
    
    double getStructBasedSeqProbRatio(MutableSequence seqI, MutableSequence seqJ, int differPosition) throws Exception{
    	double structBasedSeqProbRatio = 0;
    	
    	int[] codonArrayI = seqI.toCodonArray();
    	int[] codonArrayJ = seqJ.toCodonArray();
    	
    	int codonPosition = differPosition / 3;
    	
    	double firstOrderRatio = Math.log(inputStructure.getFirstOrderProb(codonPosition, codonArrayJ[codonPosition])) - Math.log(inputStructure.getFirstOrderProb(codonPosition, codonArrayI[codonPosition]));
    	
    	double interactionRatio = 0;
    	
		for (int m = 0; m < codonPosition; m++) {
			interactionRatio += Math.log(inputStructure.getInteractionProb(m, codonPosition, codonArrayJ[m], codonArrayJ[codonPosition])) - Math.log(inputStructure.getInteractionProb(m, codonPosition, codonArrayI[m], codonArrayI[codonPosition]));
		}
		for (int n = codonPosition + 1 ; n < codonArrayI.length; n++) {
			interactionRatio += Math.log(inputStructure.getInteractionProb(codonPosition, n, codonArrayJ[codonPosition], codonArrayJ[n])) - Math.log(inputStructure.getInteractionProb(codonPosition, n, codonArrayI[codonPosition], codonArrayI[n]));
		}
    	
    	structBasedSeqProbRatio = firstOrderRatio + interactionRatio;
    	return structBasedSeqProbRatio;
    }
    
}