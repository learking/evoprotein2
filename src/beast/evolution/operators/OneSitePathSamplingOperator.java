/**
 * 
 */
package beast.evolution.operators;

import java.util.ArrayList;
import java.util.List;

import evoprotein.evolution.datatype.MutableSequence;
import evoprotein.evolution.substitution.SubstitutionEvent;
import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.InstantHKY;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.PathBranch;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.SeqPath;
import beast.util.Randomizer;

/**
 * @author kwang2
 *
 */
@Description("One site version of pathSamplingOperator, used in every MCMC round except initialization stage")
public class OneSitePathSamplingOperator extends Operator{

	public Input<PathTree> m_pathTree = new Input<PathTree>("pathtree", "pathtree on which this operation is performed", Validate.REQUIRED);
    public Input<SiteModel.Base> m_pSiteModel = new Input<SiteModel.Base>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);
    
    int randomSite;
    int seqLength;
    int rootNr;
    int sudoRootNr;

    List<Integer> internalNodesNr;
    
    SubstitutionModel.Base m_substitutionModel;
    double fHastingsRatio;
    double [] oldSitesProb;
    double newSiteProb;

    // lambda parameter for nextExponential()
    final double lambda = 1.0;
    
    boolean firstTimeFlag = true;
    
    @Override
    public void accept(){
    	//System.out.println("Hasting ratio:" + fHastingsRatio);
    	super.accept();
    	oldSitesProb[randomSite] = newSiteProb;
    }
    
    @Override
    public void initAndValidate() throws Exception {
    	m_substitutionModel = m_pSiteModel.get().substModelInput.get();
    	
    	rootNr = m_pathTree.get().getRoot().getNr();
    	sudoRootNr = m_pathTree.get().getSudoRootNr();
    	internalNodesNr = getInternalNodesNr();
    	
    	seqLength = m_pathTree.get().getSequences().get(0).getSequence().length;
    	oldSitesProb = new double[seqLength];
    }
    
    List<Integer> getInternalNodesNr(){
    	List<Integer> internalNodes = new ArrayList<Integer>();
		for (int internalNodeIndex = m_pathTree.get().getLeafNodeCount(); internalNodeIndex < m_pathTree.get().getNodeCount(); internalNodeIndex++) {
			internalNodes.add(internalNodeIndex);
		}
		// remove the root number from internalNodeNum
		internalNodes.remove(internalNodes.indexOf(rootNr));
    	return internalNodes;
    }
    
	@Override
	public double proposal() {
		//System.out.println("#####################################################################################");

		// reset fHastingsRatio
		fHastingsRatio = 0;
		
		// register this operator with input PathTree
		PathTree pathTree = m_pathTree.get(this);
		
		if(firstTimeFlag){
			firstTimeFlag = false;		
			// let MCMC accept it
			fHastingsRatio = sampleAllSites(pathTree);
		}else{
			// after first time:
			// step 1: sample path for one site only
			randomSite = (int) (Randomizer.nextDouble() * seqLength);
			//System.out.println("************************************************************************************");
			newSiteProb = sampleOneSite(pathTree, randomSite);
			//System.out.println("************************************************************************************");
			
			if(existStopCodonThisTriplet(pathTree, randomSite)){
				fHastingsRatio = Double.NEGATIVE_INFINITY;
				System.out.println("Stop codon! fHastingsRatio set to negative infinity!");
			}else{
				//System.out.println("site number:" + randomSite);
				//System.out.println("oldSite Prob:" + oldSitesProb[randomSite] );
				//System.out.println("newSite Prob:" + newSiteProb);
				fHastingsRatio = oldSitesProb[randomSite] - newSiteProb;
			}
		}
		
		// change the tree and set the tree to be dirty
		//System.out.println("proposal Hasting ratio:" + fHastingsRatio);
		
		pathTree.setSomethingIsDirty(true);
		//System.out.println("#####################################################################################");
		return fHastingsRatio;
	}

	double sampleOneSite(PathTree pathTree, int siteNr){
		double siteProb = 0;
		siteProb += PupkoOneSite(pathTree, siteNr);
		//System.out.println("====================================================================================");
		siteProb += NielsenOneSite(pathTree, siteNr);
		//System.out.println("====================================================================================");
		return siteProb;
	}
		
	double PupkoOneSite(PathTree pathTree, int siteNr){
		double pupkoLogP = 0;
		
		// if the seq state at this site of the current node "N" is "j" and its parent seq state is "i"
		// then P(N=j | N_parent = i) is the "i*4 + j"th item in the "N"th Matrix of the list
		List<double[]> pMatrices = new ArrayList<double []>();
		for (int i=0; i< pathTree.getNodeCount(); i++){
			pMatrices.add(new double[4*4]);
		}
		
		/* deal with leaf nodes first (since they can be readily calculated) */
		for (int leafNodeNr = 0; leafNodeNr < pathTree.getLeafNodeCount(); leafNodeNr++) {
			setLeafpMatrix(pMatrices.get(leafNodeNr), pathTree.getNode(leafNodeNr).getLength(), pathTree.getSequences().get(leafNodeNr).getSequence()[siteNr]);
		}
		
		// create a list to store numbers for internal nodes
		List<Integer> internalNodeIndices = new ArrayList<Integer>();
		for (int internalNodeIndex = pathTree.getLeafNodeCount(); internalNodeIndex < pathTree.getNodeCount(); internalNodeIndex++) {
			internalNodeIndices.add(internalNodeIndex);
		}
		// remove the root number from internalNodeNum
		internalNodeIndices.remove(internalNodeIndices.indexOf(rootNr));
		
		while(internalNodeIndices.size() > 0 ){
			for (int i = 0; i < internalNodeIndices.size(); i++) {
				int currentNodeNr = internalNodeIndices.get(i);
				Node currentNode = pathTree.getNode(currentNodeNr);
				// System.out.println("internalnode:" + currentNodeNr);
				if (currentNode.getLength() > 0) {
					int leftNum = currentNode.getLeft().getNr();
					int rightNum = currentNode.getRight().getNr();
				    // if pMatrices of both left and right child node has been calculated
					if ((getpMatrixSum(pMatrices.get(leftNum))!=0)&&(getpMatrixSum(pMatrices.get(rightNum))!=0)) {
						// set up pMatrix for this internal node
						setInternalpMatrix(pMatrices.get(currentNodeNr), pMatrices.get(leftNum), pMatrices.get(rightNum), currentNode.getLength());
						internalNodeIndices.remove(i);
					}
				}
				
				if (currentNode.getNr() == sudoRootNr){
					int leftNum = currentNode.getLeft().getNr();
					int rightNum = currentNode.getRight().getNr();
				    // if pMatrices of both left and right child node has been calculated
					if ((getpMatrixSum(pMatrices.get(leftNum))!=0)&&(getpMatrixSum(pMatrices.get(rightNum))!=0)) {
						// get tip node that is the immediate child of root
						int tipNodeNr = -1;
						List<Node> rootChildren = pathTree.getRoot().getChildren();
						for (Node tmpNode : rootChildren) {
							if (tmpNode.isLeaf()) {
								tipNodeNr = tmpNode.getNr();
							}
						}
						// too many parameters (to be simplified later)
						pupkoLogP += setSudoRoot(pathTree, siteNr, pMatrices.get(leftNum), pMatrices.get(rightNum), pMatrices.get(tipNodeNr));
						// trace back and set all internal nodes based on pMatrices
						pupkoLogP += setInternalNodes(pathTree, siteNr, pMatrices);
						internalNodeIndices.remove(i);
					}
				}
				
			}
		}

		// let the root have the same sequence as the sudoRoot
		MutableSequence newRootSeq = pathTree.getSequences().get(sudoRootNr).copy();
		pathTree.getSequences().get(rootNr).setSequence(newRootSeq.getSequence());
		
		return pupkoLogP;
	}
	
	public double setSudoRoot(PathTree pathTree, int seqSite, double [] leftpMatrix, double [] rightpMatrix, double[] tipNodepMatrix){
		
		double sudoRootNucleoStateProb = 0.0;
		final double [] frequences = m_substitutionModel.getFrequencies();
		double [] nucleoStateProbs = new double [4];
		for(int nucleoState=0; nucleoState < 4; nucleoState++){
			double nucleoStateFreq = frequences[nucleoState];
			double leftpMatrixRowSum = getpMatrixRowSum(leftpMatrix, nucleoState);
			double rightpMatrixRowSum = getpMatrixRowSum(rightpMatrix, nucleoState);
			double tipNodepMatrixRowSum = getpMatrixRowSum(tipNodepMatrix, nucleoState);
			nucleoStateProbs[nucleoState] = nucleoStateFreq * leftpMatrixRowSum * rightpMatrixRowSum * tipNodepMatrixRowSum;
		}
		
		// normalize nucleoStateFreq using its added-up total
		double totalNucleoStateProb = 0;
		for (int i=0; i < nucleoStateProbs.length; i++) {
			totalNucleoStateProb += nucleoStateProbs[i];
		}

		// compute CDF based on nucleoStateProbs
		double [] nucleoStateCDF = new double [4];
		double cumulativeProb = 0;
		for (int i=0; i < nucleoStateProbs.length; i++) {
			nucleoStateProbs[i] = nucleoStateProbs[i] / totalNucleoStateProb;
			cumulativeProb += nucleoStateProbs[i];
			nucleoStateCDF[i] = cumulativeProb;
		}
		
		// pick a nucleoState randomly based on nucleoStateProbs
		int sudoRootNucleoState = Randomizer.randomChoice(nucleoStateCDF);
		sudoRootNucleoStateProb = nucleoStateProbs[sudoRootNucleoState];
		pathTree.getSequences().get(sudoRootNr).getSequence()[seqSite] = sudoRootNucleoState;
		
		return Math.log(sudoRootNucleoStateProb);
	}
	
	public double setInternalNodes(PathTree pathTree, int seqSite, List<double []> pMatrices) {
		double internalNodesLogP = 0;
		Node leftNode = pathTree.getNode(sudoRootNr).getLeft();
		Node rightNode = pathTree.getNode(sudoRootNr).getRight();
		// traverseInternalNode is a recursive function
		traverseInternalNode(internalNodesLogP, leftNode, pathTree, seqSite, pMatrices);
		traverseInternalNode(internalNodesLogP, rightNode, pathTree, seqSite, pMatrices);
		return internalNodesLogP;
	}
	
	public void traverseInternalNode(double internalNodesLogP, Node node, PathTree pathTree, int seqSite, List<double []> pMatrices) {
		if(!node.isLeaf()) {
			double [] pMatrix = pMatrices.get(node.getNr());
			int parentNucleoState = pathTree.getSequences().get(node.getParent().getNr()).getSequence()[seqSite];
			
			double [] currentpVector = getpMatrixOneRow(pMatrix, parentNucleoState);
			double currentpVectorSum = getpMatrixRowSum(pMatrix, parentNucleoState);
			for (int i = 0; i < currentpVector.length ; i++) {
				currentpVector[i] = currentpVector[i] / currentpVectorSum;
			}
			
			// compute CDF based on nucleoStateProbs
			double [] thisNodeStateCDF = new double [4];
			double cumulativeProb = 0;
			for (int i=0; i < currentpVector.length; i++) {
				cumulativeProb += currentpVector[i];
				thisNodeStateCDF[i] = cumulativeProb;
			}
			
			// pick a nucleoState randomly based on nucleoStateProbs
			int thisNodeState = Randomizer.randomChoice(thisNodeStateCDF);
			double thisNodeStateProb = currentpVector[thisNodeState];
			internalNodesLogP += Math.log(thisNodeStateProb);
			pathTree.getSequences().get(node.getNr()).getSequence()[seqSite] = thisNodeState;
			
			//recursive part:
			traverseInternalNode(internalNodesLogP, node.getLeft() , pathTree, seqSite, pMatrices);
			traverseInternalNode(internalNodesLogP, node.getRight() , pathTree, seqSite, pMatrices);
			
		}
	}
	
	public void setLeafpMatrix(double [] pMatrix, double branchLength, int nucleotideState){
		double [] transitionMatrix = new double [4*4];
		m_substitutionModel.getTransitionProbabilities(null, branchLength, 0, 1.0, transitionMatrix);
		for(int i=0; i<4; i++){
			pMatrix[i*4 + nucleotideState] = transitionMatrix[i*4 + nucleotideState];
		}
	}

	public void setInternalpMatrix(double [] currentpMatrix, double [] leftpMatrix, double [] rightpMatrix, double branchLength){
		double [] transitionMatrix = new double [4*4];
		m_substitutionModel.getTransitionProbabilities(null, branchLength, 0, 1.0, transitionMatrix);

		for(int parentNodeNucleoState = 0; parentNodeNucleoState < 4; parentNodeNucleoState++){
			for(int currentNodeNucleoState = 0; currentNodeNucleoState < 4; currentNodeNucleoState++){
				double transitionProb = transitionMatrix[parentNodeNucleoState * 4 + currentNodeNucleoState];
				double leftpMatrixRowSum = getpMatrixRowSum(leftpMatrix, currentNodeNucleoState);
				double rightpMatrixRowSum = getpMatrixRowSum(rightpMatrix, currentNodeNucleoState);
				currentpMatrix[parentNodeNucleoState * 4 + currentNodeNucleoState] = transitionProb * leftpMatrixRowSum * rightpMatrixRowSum;
			}
		}
	}
	
	public double getpMatrixRowSum(double [] pMatrix, int rowNr){

		double rowSum = 0.0;
		int startPosition = rowNr * 4;
		int endPosition = startPosition + 4;
		for (int i = startPosition; i< endPosition; i++){
			rowSum += pMatrix[i];
		}

		return rowSum;
	}
	
	public double getpMatrixSum(double [] pMatrix) {
		double sum = 0;
		for (double p : pMatrix) {
			sum += p;
		}
		return sum;
	}
	
	public double [] getpMatrixOneRow(double [] pMatrix, int rowNr) {
		double [] thisRow = new double [4];
		int thisRowNr = rowNr;
		int startPosition = thisRowNr * 4;
		int endPosition = startPosition + 4;
		for (int i = startPosition; i < endPosition; i++) {
			int thisRowIndex = i - startPosition;
			thisRow[thisRowIndex] = pMatrix[i];
		}
		return thisRow;
	}
	
	double NielsenOneSite(PathTree pathTree, int siteNr){
		double nielsenLogP = 0;
		
		for (int branchNr = 0; branchNr < pathTree.getBranches().size(); branchNr++) {
			if((branchNr != rootNr) && (branchNr != sudoRootNr)) {
				try {
					PathBranch thisBranch = pathTree.getBranch(branchNr);
					double thisBranchLength = pathTree.getNode(thisBranch.getEndNodeNr()).getLength();
					int endNodeNr = thisBranch.getEndNodeNr();
					int beginNodeNr = thisBranch.getBeginNodeNr();
					MutableSequence childSeq = pathTree.getSequences().get(endNodeNr);
					MutableSequence parentSeq = pathTree.getSequences().get(beginNodeNr);
					int parentNucleoState = parentSeq.getSequence()[siteNr];
					int childNucleoState = childSeq.getSequence()[siteNr];
					 //System.out.println("branch nr:" + branchNr + " parent nucleo:" + parentNucleoState + " child nucleo:" + childNucleoState);
					nielsenLogP += NielsenSampleOneBranchOneSite(thisBranch, siteNr, thisBranchLength, parentNucleoState, childNucleoState);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		return nielsenLogP;
	}
	
	public double NielsenSampleOneBranchOneSite(PathBranch thisBranch, int seqSite, double thisBranchLength, int parentNucleoState, int childNucleoState){
		double nielsenOneBranchOneSiteLogP = 0;
		
		final int childState = childNucleoState;
		final int parentState = parentNucleoState;
		final double totalTime = thisBranchLength;

		List<SubstitutionEvent> substitutionEvents = new ArrayList<SubstitutionEvent>();
		
		// CDFs for sampling another different nucleotide
		final double [][] differentCDFs = getDifferentNucleoCDF();
		
		// parameters used in the calculation
		int lastNucleotide = -1;
		double currentTime = 0;
		int currentState = parentState;
		int beginNucleotide = -1;
		int endNucleotide = -1;
		double timeInterval = 0;
		
		if(parentState == childState){
			while(lastNucleotide != childState){
				// initialize stuff
				substitutionEvents.clear();			
				currentTime = 0;
				currentState = parentState;
				
				while(currentTime < totalTime){
					timeInterval = Randomizer.nextExponential(lambda);
					if(currentTime + timeInterval >= totalTime){
						currentTime = totalTime;
						lastNucleotide = currentState;
					}else{
						currentTime += timeInterval;
						beginNucleotide = currentState;
						currentState = Randomizer.randomChoice(differentCDFs[currentState]);
						endNucleotide = currentState;
						// create new substitutionevent and add it to events
						// System.out.println("begin:" + beginNucleotide + " end:" + endNucleotide);
						substitutionEvents.add(new SubstitutionEvent(beginNucleotide, endNucleotide, timeInterval));
					}
				}
			}

			//System.out.println("# of subs:" + substitutionEvents.size());
			// copy substitutionEvents to ...
			if(substitutionEvents.size() != 0){
				thisBranch.setMutationPath(seqSite, substitutionEvents);
			}else{
				thisBranch.getMutationPath(seqSite).clear();
			}
			
		}
		else{
			while (lastNucleotide != childState) {
				// initialize stuff
				substitutionEvents.clear();			
				currentTime = 0;
				currentState = parentState;
				
				boolean firstSample = true;
				while (currentTime < totalTime) {
					if (firstSample) {
						timeInterval = sampleFirstSubstitutionTime(totalTime);
						currentTime = currentTime + timeInterval;
						
						beginNucleotide = currentState;
						currentState = Randomizer.randomChoice(differentCDFs[currentState]);
						endNucleotide = currentState;
						substitutionEvents.add(new SubstitutionEvent(beginNucleotide, endNucleotide, timeInterval));
						
						firstSample = false;
					} else {
						timeInterval = Randomizer.nextExponential(lambda);
						if (currentTime + timeInterval >= totalTime) {
							currentTime = totalTime;
							lastNucleotide = currentState;
						} else {
							currentTime += timeInterval;
							
							beginNucleotide = currentState;
							currentState = Randomizer.randomChoice(differentCDFs[currentState]);
							endNucleotide = currentState;
							substitutionEvents.add(new SubstitutionEvent(beginNucleotide, endNucleotide, timeInterval));
						}
					}
				}
			}
			
			// copy substitutionEvents to ...
			if(substitutionEvents.size() != 0){
				thisBranch.setMutationPath(seqSite, substitutionEvents);
			}else{
				System.err.print("Error: expecting at least one substitution at this site on this branch!");
			}
		}

		nielsenOneBranchOneSiteLogP += calculateOneBranchOneSiteLogP(childState, parentState, totalTime, substitutionEvents);

		return nielsenOneBranchOneSiteLogP;
	}
	
	public double [][] getDifferentNucleoCDF(){
		final double [] frequences = m_substitutionModel.getFrequencies();
		double [][] differentNucleoCDFs = new double [4][4];
		
		double [] CDFwithoutA = getCDFwithoutOneNucleo(frequences, 0);
		double [] CDFwithoutC = getCDFwithoutOneNucleo(frequences, 1);
		double [] CDFwithoutG = getCDFwithoutOneNucleo(frequences, 2);
		double [] CDFwithoutT = getCDFwithoutOneNucleo(frequences, 3);
		
		differentNucleoCDFs[0]= CDFwithoutA;
		differentNucleoCDFs[1]= CDFwithoutC;
		differentNucleoCDFs[2]= CDFwithoutG;
		differentNucleoCDFs[3]= CDFwithoutT;
		
		return differentNucleoCDFs;
	}
	
	public double [] getCDFwithoutOneNucleo(double [] frequences, int nucleoState){
		double [] tmpFreqs = new double[4];
		double [] tmpCDF = new double [4];
		
		for (int i = 0; i < 4 ; i++){
			if(i == nucleoState){
				tmpFreqs[i] = 0;
			}else{
				tmpFreqs[i] = frequences[i];
			}
		}
		
		double [] tmpPDF = Randomizer.getNormalized(tmpFreqs);
		double cumulativeProb = 0;
		for (int i = 0 ; i < 4 ; i++) {
			cumulativeProb += tmpPDF[i];
			tmpCDF[i] = cumulativeProb;
		}
		return tmpCDF;
	}
	
	public double calculateOneBranchOneSiteLogP(int childState, int parentState, double branchLength, List<SubstitutionEvent> substitutionEvents) {

		double oneBranchOneSiteLogP = 0.0;
		double[] instantMatrix = new double[4 * 4];
		// this needs to be changed to more general method
		((InstantHKY) m_substitutionModel).getInstantRateMatrix(instantMatrix);

		// get "end" and "begin" nucleotide
		int endNucleotide = childState;
		int beginNucleotide = parentState;

		int transitionOutEventCode;
		int transitionEventCode;
		double currentBranchLength = branchLength;
		List<SubstitutionEvent> currentSubstitutionEvents = substitutionEvents;

		if (currentSubstitutionEvents.size() == 0) {
			// diagonal elements are already negative
			transitionOutEventCode = beginNucleotide * 4 + beginNucleotide;
			oneBranchOneSiteLogP += instantMatrix[transitionOutEventCode]
					* currentBranchLength;
		} else {
			// ...
			double cumulativeHeight = 0.0;
			int lastNucleotide = -1;
			for (int substitutionIndex = 0; substitutionIndex < currentSubstitutionEvents
					.size(); substitutionIndex++) {

				int previousNucleotide = currentSubstitutionEvents.get(
						substitutionIndex).getPreviousNucleotide();
				int currentNucleotide = currentSubstitutionEvents.get(
						substitutionIndex).getCurrentNucleotide();
				// determine the last nucleotide
				if (substitutionIndex == (currentSubstitutionEvents.size() - 1)) {
					lastNucleotide = currentNucleotide;
				}

				double currentTimeInterval = currentSubstitutionEvents.get(
						substitutionIndex).getTimeInterval();
				cumulativeHeight += currentTimeInterval;
				transitionOutEventCode = previousNucleotide * 4
						+ previousNucleotide;
				transitionEventCode = previousNucleotide * 4
						+ currentNucleotide;
				oneBranchOneSiteLogP += Math.log(instantMatrix[transitionEventCode])
						+ instantMatrix[transitionOutEventCode]
						* currentTimeInterval;
			}

			double unchangedHeight = currentBranchLength - cumulativeHeight;

			transitionOutEventCode = lastNucleotide * 4 + lastNucleotide;
			oneBranchOneSiteLogP += instantMatrix[transitionOutEventCode]
					* unchangedHeight;
		}

		double [] transitionMatrix = new double [4*4];
		m_substitutionModel.getTransitionProbabilities(null, currentBranchLength, 0, 1, transitionMatrix);
		oneBranchOneSiteLogP = oneBranchOneSiteLogP - Math.log(transitionMatrix[beginNucleotide*4 + endNucleotide]);
		return oneBranchOneSiteLogP;
	}
	
	public double sampleFirstSubstitutionTime(double totalTime) {
		// store the final result
		double sampledFirstTime;

		// draw a random number between 0 and 1
		double randNum = Randomizer.nextDouble();

		// calculate the time for the first substitution
		sampledFirstTime = -Math
				.log(1.0 - randNum * (1 - Math.exp(-totalTime)));
		return sampledFirstTime;
	}
	
	
	
	double sampleAllSites(PathTree pathTree){
		
		for (int startSite = 0; startSite < (seqLength - 2) ; startSite = startSite + 3) {
				oldSitesProb[startSite] = sampleOneSite(pathTree, startSite);
				oldSitesProb[startSite + 1] = sampleOneSite(pathTree, startSite + 1);
				oldSitesProb[startSite + 2] = sampleOneSite(pathTree, startSite + 2);
			while (existStopCodonThisTriplet(pathTree, startSite)) {
				oldSitesProb[startSite] = sampleOneSite(pathTree, startSite);
				oldSitesProb[startSite + 1] = sampleOneSite(pathTree, startSite + 1);
				oldSitesProb[startSite + 2] = sampleOneSite(pathTree, startSite + 2);
			}
			
		}
		
		return Double.POSITIVE_INFINITY;
	}
	
	// checkers
	boolean existStopCodonThisTriplet(PathTree pathTree, int siteNr){
		boolean stopCodonFlag = false;
		int startSite = siteNr - ( siteNr % 3 );
		
		// check triplet type for internal nodes
		for (Integer internalNodeIndex : internalNodesNr) {
			MutableSequence currentSeq = pathTree.getSequences().get(internalNodeIndex.intValue());
			if(currentSeq.getNucleotide(startSite) == 3){
				if(currentSeq.getNucleotide(startSite + 1) == 0){
					if(currentSeq.getNucleotide(startSite + 2) == 0){
						stopCodonFlag = true;
						break;						
					}
					if(currentSeq.getNucleotide(startSite + 2) == 2){
						stopCodonFlag = true;
						break;						
					}					
				}
				if((currentSeq.getNucleotide(startSite + 1) == 2) && (currentSeq.getNucleotide(startSite + 2) == 0)){
					stopCodonFlag = true;
					break;						
				}
			}
		}
		
		// if no stop codon within internal nodes
		if(!stopCodonFlag){
		// check triplet type for seq path along each branch
			for (int branchNr = 0; branchNr < pathTree.getBranches().size(); branchNr++) {
				if((branchNr != rootNr) && (branchNr != sudoRootNr)) {
					try {
						PathBranch thisBranch = pathTree.getBranch(branchNr);
						int endNodeNr = thisBranch.getEndNodeNr();
						int beginNodeNr = thisBranch.getBeginNodeNr();
						MutableSequence childSeq = pathTree.getSequences().get(endNodeNr);
						MutableSequence parentSeq = pathTree.getSequences().get(beginNodeNr);
						SeqPath codonSeqPath = thisBranch.getCodonSeqPath(siteNr, parentSeq, childSeq);
						if(codonSeqPath.existStopCodon()){
							stopCodonFlag = true;
							break;
						}
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		
		return stopCodonFlag;
	}
	
	
	
}
