package beast.evolution.substitutionmodel;

import evoprotein.evolution.datatype.CodonUtil;
import evoprotein.evolution.datatype.MutableSequence;
import evoprotein.proteinstructure.InputStructure;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.GappedAlignment;
import beast.evolution.datatype.DataType;
import beast.evolution.tree.Node;

public class ProteinCodingDNASubstModel extends CalculationNode {

	// determine inputs	
	public Input<RealParameter> kappa = new Input<RealParameter>("kappa", "kappa parameter in HKY model", Validate.REQUIRED);
	
	// for use in neutral model
	public Input<Frequencies> frequenciesInput =
            new Input<Frequencies>("frequencies", "our substitution model equilibrium state frequencies", Validate.REQUIRED);
	// for use in struct based model
	public Input<InputStructure> m_inputStructure = new Input<InputStructure>("inputStructure", "input Structure pre-calculated", Validate.REQUIRED);

	// add gapped alignment as input so that we would be able to know where deletion-caused gaps happen
	public Input<Alignment> m_alignment = new Input<Alignment>("alignment", "gapped alignment", Validate.REQUIRED);
	
	static CodonUtil codonUtil = new CodonUtil();
	
	// define variables here
    Frequencies frequencies;
    InputStructure inputStructure;
    double freqA;
    double freqC;
    double freqG;
    double freqT;
    double Y;
    double scalingFactor;
    int interactionRange;
    
	// initAndValidate
    @Override
    public void initAndValidate() throws Exception{
    	frequencies = frequenciesInput.get();
    	inputStructure = m_inputStructure.get();
    	
    	// cast to gapped alignment
    	GappedAlignment alignment = (GappedAlignment) m_alignment.get();
    	
    	// remove terms associated with deletion-caused gaps
    	inputStructure.removeGapRelatedTerms(alignment.getDeletionPositions());
    	
    	interactionRange = 10;
    }
    
    public double getSubstAwayRate(MutableSequence seqI) throws Exception{
		//System.out.println("Rate away start:" + System.currentTimeMillis());
    	double awayRate = 0;
    	int[] codonArrayI = seqI.toCodonArray();
    	
    	//MutableSequence tmpSeq;
    	int[] originalNucleoSeq = seqI.getSequence();
    	MutableSequence tmpSeq = new MutableSequence(originalNucleoSeq.length);
    	
    	for (int site = 0 ; site < seqI.getSequence().length; site++) {
    		// init tmpSeq when dealing with each site
    		//tmpSeq = seqI.copy();
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

		// find where these two sequences differ (both location and value)
		//int differPosition = getDifferPosition(seqI, seqJ);

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
    	//int differCodon = getDifferCodon(seqJ, codonDifferPosition);
    	double structBasedSeqProbRatio = getStructBasedSeqProbRatio(codonArrayI, differCodon, codonDifferPosition);

    	logTAU = structBasedSeqProbRatio - neutralSeqProbRatio;
    	return logTAU;
    }
    
    int getDifferCodon(MutableSequence seqJ, int codonDifferPosition){
    	int differCodon = -1;
    	differCodon = codonUtil.translate(seqJ, codonDifferPosition * 3);
    	return differCodon;
    }
    
    // full interaction version of getStructBasedSeqProbRatio
    /*
    double getStructBasedSeqProbRatio(int[] codonArrayI, int differCodon, int codonDifferPosition) throws Exception{
    	double structBasedSeqProbRatio = 0;
    	
    	double firstOrderRatio = inputStructure.getFirstOrderLogProb(codonDifferPosition, differCodon) - inputStructure.getFirstOrderLogProb(codonDifferPosition, codonArrayI[codonDifferPosition]);
    	
    	double interactionRatio = 0;
    	
    	// here, m refers to a codon position, it should in
    	
		for (int m = 0; m < codonDifferPosition; m++) {
			interactionRatio += inputStructure.getInteractionLogProb(m, codonDifferPosition, codonArrayI[m], differCodon) - inputStructure.getInteractionLogProb(m, codonDifferPosition, codonArrayI[m], codonArrayI[codonDifferPosition]);
		}
		
		for (int n = codonDifferPosition + 1 ; n < codonArrayI.length; n++) {
			interactionRatio += inputStructure.getInteractionLogProb(codonDifferPosition, n, differCodon, codonArrayI[n]) - inputStructure.getInteractionLogProb(codonDifferPosition, n, codonArrayI[codonDifferPosition], codonArrayI[n]);
		}	
    	
    	structBasedSeqProbRatio = firstOrderRatio + interactionRatio;
    	//System.out.println(interactionRatio);
    	return structBasedSeqProbRatio;
    }
    */
    
    // range 10 version of getStructBasedSeqProbRatio
    double getStructBasedSeqProbRatio(int[] codonArrayI, int differCodon, int codonDifferPosition) throws Exception{
    	double structBasedSeqProbRatio = 0;
    	
    	double firstOrderRatio = inputStructure.getFirstOrderLogProb(codonDifferPosition, differCodon) - inputStructure.getFirstOrderLogProb(codonDifferPosition, codonArrayI[codonDifferPosition]);
    	//System.out.println("first order:" + firstOrderRatio);
    	double interactionRatio = 0;
    	
    	// here, m refers to a codon position, it should in
    	
		for (int m = getInteractionRangeLeftBound(codonDifferPosition); m < codonDifferPosition; m++) {
			interactionRatio += inputStructure.getInteractionLogProb(m, codonDifferPosition, codonArrayI[m], differCodon) - inputStructure.getInteractionLogProb(m, codonDifferPosition, codonArrayI[m], codonArrayI[codonDifferPosition]);
		}
		
		for (int n = codonDifferPosition + 1 ; n < getInteractionRangeRightBound(codonDifferPosition, codonArrayI.length); n++) {
			interactionRatio += inputStructure.getInteractionLogProb(codonDifferPosition, n, differCodon, codonArrayI[n]) - inputStructure.getInteractionLogProb(codonDifferPosition, n, codonArrayI[codonDifferPosition], codonArrayI[n]);
		}	
    	
    	structBasedSeqProbRatio = firstOrderRatio + interactionRatio;
    	//System.out.println(interactionRatio);
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
    
    public double getRootSeqLogP(int[] rootCodonSeq) {
    	double rootCodonSeqLogP = 0.0;
    	
    	// deal with first order terms
    	for(int codonPosition = 0; codonPosition < rootCodonSeq.length; codonPosition++) {
    		rootCodonSeqLogP += inputStructure.getFirstOrderLogProb(codonPosition, rootCodonSeq[codonPosition]);
    	}
    	
    	// deal with second order terms
    	for(int firstCodonPosition = 0; firstCodonPosition < (rootCodonSeq.length-1); firstCodonPosition++) {
    		for(int secondCodonPosition = (firstCodonPosition+1); secondCodonPosition < rootCodonSeq.length; secondCodonPosition++) {
    			rootCodonSeqLogP += inputStructure.getInteractionLogProb(firstCodonPosition, secondCodonPosition, rootCodonSeq[firstCodonPosition], rootCodonSeq[secondCodonPosition]);
    		}   		
    	}
    	
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
