/**
 * 
 */
package beast.evolution.substitutionmodel;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.GappedAlignment;
import evoprotein.evolution.datatype.CodonUtil;
import evoprotein.evolution.datatype.MutableSequence;
import evoprotein.proteinstructure.InputStructure;
import evoprotein.proteinstructure.InputTwoStruct;

/**
 * @author kuangyu
 *
 */
public class TwoStructSubstModel extends CalculationNode {
	
	// determine inputs	
	public Input<RealParameter> kappa = new Input<RealParameter>("kappa", "kappa parameter in HKY model", Validate.REQUIRED);
	
	// newly introduced "f" parameter
	public Input<RealParameter> f = new Input<RealParameter>("f", "f parameter", Validate.REQUIRED);	
	
	// for use in neutral model
	public Input<Frequencies> m_frequencies =
            new Input<Frequencies>("frequencies", "our substitution model equilibrium state frequencies", Validate.REQUIRED);
	// for use in struct based model
	// this old version of inputStructure can only handle one structure case
	//public Input<InputStructure> m_inputStructure = new Input<InputStructure>("inputStructure", "input Structure pre-calculated", Validate.REQUIRED);
	public Input<InputTwoStruct> m_inputTwoStruct = new Input<InputTwoStruct>("inputTwoStruct", "two input Structures pre-calculated", Validate.REQUIRED);
	
	// add gapped alignment as input so that we would be able to know where deletion-caused gaps happen
	public Input<Alignment> m_alignment = new Input<Alignment>("alignment", "gapped alignment", Validate.REQUIRED);
	
	static CodonUtil codonUtil = new CodonUtil();
	
	// define variables here
    Frequencies frequencies;
    // differ from one struct
    InputTwoStruct inputTwoStruct;
    double freqA;
    double freqC;
    double freqG;
    double freqT;
    double Y;
    double scalingFactor;
    int interactionRange;

    // initAndValidate
    @Override
    public void initAndValidate(){
    	frequencies = m_frequencies.get();
    	inputTwoStruct = m_inputTwoStruct.get();
    	
    	// cast to gapped alignment
    	GappedAlignment alignment = (GappedAlignment) m_alignment.get();
    	// do sth here to handle deletion-caused gaps
    	inputTwoStruct.removeGapRelatedTerms(alignment.getDeletionPositions());
    	
    	// range is set here
    	interactionRange = 10;
    }
    
    public double getSubstAwayRate(MutableSequence seqI) throws Exception{
		//System.out.println("Rate away start:" + System.currentTimeMillis());
    	double awayRate = 0;
    	int[] originalNucleoSeq = seqI.getSequence();
    	int seqLength = originalNucleoSeq.length;
    	int[] codonArrayI = seqI.toCodonArray();
    	
    	//MutableSequence tmpSeq;
    	MutableSequence tmpSeq = new MutableSequence(originalNucleoSeq.length);
    	
    	for (int site = 0 ; site < seqLength; site++) {
    		// init tmpSeq when dealing with each site
    		tmpSeq.setSequence(originalNucleoSeq);
    		int originalNucleotide = seqI.getSequence()[site];
    		int[] changableNucleotides = getChangableNucleotides(originalNucleotide);
    		for (int i = 0; i < changableNucleotides.length; i++) {
    			tmpSeq.mutate(site, changableNucleotides[i]);
    			//System.out.println(tmpSeq.toString());
    			int mutatedCodonStartSite = site - (site%3);
    			int differCodon = codonUtil.translate(tmpSeq, mutatedCodonStartSite);
    			if(differCodon != -1) {
    				awayRate += getSubstitutionRate(seqI, tmpSeq, codonArrayI, site, differCodon);
    			}
    		}
    	}
		//System.out.println("Rate away end:" + System.currentTimeMillis());
    	return awayRate;
    }    

    
    int[] getChangableNucleotides(int originalNucleotide){
    	int [] changableNucleotides;
    	switch(originalNucleotide){
    	case 0:
    		changableNucleotides = new int[]{1,2,3};
    		break;
    	case 1:
    		changableNucleotides = new int[]{0,2,3};
    		break;
    	case 2:
    		changableNucleotides = new int[]{0,1,3};
    		break;
    	case 3:
    		changableNucleotides = new int[]{0,1,2};
    		break;
    	default:
    		changableNucleotides = new int[]{-1,-1,-1};
    	} 
    	return changableNucleotides;
    }
    
	public double getSubstitutionRate(MutableSequence seqI, MutableSequence seqJ, int[] codonArrayI, int differPosition, int differCodon) throws Exception{
		
		scalingFactor = getScalingFactor();
				
		double substitutionRate = 0;

		double logTAU = getLogTAU(seqI, seqJ, codonArrayI, differPosition, differCodon);
		
		substitutionRate = logTAU / (1 - 1/Math.exp(logTAU));
		
		if(isTransition(seqI.getNucleotide(differPosition), seqJ.getNucleotide(differPosition))){
			double k = kappa.get().getValue();
			substitutionRate = scalingFactor * substitutionRate * frequencies.getFreqs()[seqJ.getNucleotide(differPosition)] * k;
		}else{
			substitutionRate = scalingFactor * substitutionRate * frequencies.getFreqs()[seqJ.getNucleotide(differPosition)];
		}

		return substitutionRate;
	}
    
	double getScalingFactor(){
		double u;
		double freqSquareSum = 0;
		double[] freqs = frequencies.getFreqs();
		for (int i = 0 ; i < freqs.length ; i++) {
			freqSquareSum += freqs[i] * freqs[i];
		}
		u = 1.0 / (1.0 - freqSquareSum);
		return u; 
	}
    
    boolean isTransition(int firstNucleotide, int secondNucleotide){
    	int sum = firstNucleotide + secondNucleotide;
    	if(sum == 2 || sum == 4){
    		return true;
    	}else{
    		return false;
    	}
    }
    
    double getLogTAU(MutableSequence seqI, MutableSequence seqJ, int[] codonArrayI, int differPosition, int differCodon) throws Exception{
    	double logTAU;
    	double neutralSeqProbRatio = Math.log(frequencies.getFreqs()[seqJ.getSequence()[differPosition]] / frequencies.getFreqs()[seqI.getSequence()[differPosition]]);

    	int codonDifferPosition = differPosition / 3;
    	
    	// if non-Synonymous substitution
    	double structBasedSeqProbRatio;
    	if(codonArrayI[codonDifferPosition] != differCodon ){
    		structBasedSeqProbRatio = getStructBasedSeqProbRatio(codonArrayI, differCodon, codonDifferPosition);
    	}else{
    		// if synonymous substitution
    		structBasedSeqProbRatio = 0;
    	}

    	logTAU = structBasedSeqProbRatio - neutralSeqProbRatio;
    	return logTAU;
    }
    
    int getDifferCodon(MutableSequence seqJ, int codonDifferPosition){
    	int differCodon = -1;
    	differCodon = codonUtil.translate(seqJ, codonDifferPosition * 3);
    	return differCodon;
    }
    
    // range 10 version of getStructBasedSeqProbRatio
    // change it to allow inputStructure as input
    double getStructBasedSeqProbRatio(int[] codonArrayI, int differCodon, int codonDifferPosition) throws Exception{
    	double fNow = f.get().getValue();
    	double structBasedSeqProbRatio = 0;
    	
    	double firstOrderRatio = 0;
    	// get firstOrderRatio for both structures, f*firstOrderStructA + (1-f)*firstOrderStructB
    	firstOrderRatio = inputTwoStruct.getFirstOrderRatio(codonDifferPosition, differCodon, codonArrayI[codonDifferPosition], fNow);

    	// our new way to calculate interactionRatio
    	double interactionRatio = 0;
    	int leftBound = getInteractionRangeLeftBound(codonDifferPosition);
    	int rightBound = getInteractionRangeRightBound(codonDifferPosition, codonArrayI.length);
    	interactionRatio = inputTwoStruct.getInteractionRatio(codonArrayI, leftBound, rightBound, codonDifferPosition, differCodon, fNow);
    	
    	structBasedSeqProbRatio = firstOrderRatio + interactionRatio;
    	return structBasedSeqProbRatio;
    }
    
    int getInteractionRangeLeftBound(int codonDifferPosition){
    	int leftBound = -1;
    	if (codonDifferPosition <= interactionRange) {
    		leftBound = 0;
    	}else{
    		leftBound = codonDifferPosition - interactionRange;
    	}
    	return leftBound;
    }

    int getInteractionRangeRightBound(int codonDifferPosition, int codonArrayLength){
    	int rightBound = -1;
    	if (codonDifferPosition + interactionRange + 1 >= codonArrayLength) {
    		rightBound = codonArrayLength;
    	}else{
    		rightBound = codonDifferPosition + interactionRange + 1;
    	}
    	return rightBound;
    }
    
    // ????? find original expression for rootSeqLogProb and then decide how to deal with it
    public double getRootSeqLogP(int[] rootCodonSeq) {
    	double fNow = f.get().getValue();
    	double rootCodonSeqLogP = 0.0;
    	
    	rootCodonSeqLogP = inputTwoStruct.getRootSeqLogP(rootCodonSeq, fNow);
    	
    	return rootCodonSeqLogP;
    }
    
    
    /**
     * CalculationNode implementations *
     */
    @Override
    protected boolean requiresRecalculation() {
    	// is this statement true???????????????
        // we only get here if something is dirty
        return true;
    }
	
}
