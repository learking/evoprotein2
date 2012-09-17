/**
 * 
 */
package beast.evolution.operators;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import evoprotein.evolution.datatype.MutableSequence;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.PathTree;
import beast.util.Randomizer;

/**
 * @author kuangyu
 *
 */


public class PathSamplingOperator extends Operator {
	
	public Input<PathTree> m_pathTree = new Input<PathTree>("pathtree", "pathtree on which this operation is performed", Validate.REQUIRED);
    public Input<SiteModel.Base> m_pSiteModel = new Input<SiteModel.Base>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);
	
    SubstitutionModel.Base m_substitutionModel;
	
    @Override
    public void initAndValidate() throws Exception {
    	m_substitutionModel = m_pSiteModel.get().m_pSubstModel.get();
    }
    
	/*
	 * @see beast.core.Operator#proposal()
	 */
	@Override
	public double proposal() {
		
		// register this operator with input PathTree
		PathTree pathTree = m_pathTree.get(this);
		int rootNr = pathTree.getRoot().getNr();
		int sudoRootNr = 0;
		for (Node childNode : pathTree.getRoot().getChildren()) {
			if(!childNode.isLeaf()) {
				sudoRootNr = childNode.getNr();
			}
		}

		// Is there any way to utilize old path density without recalculation?
		// for example, if the proposal gets rejected
		double oldPathDensity, newPathDensity, fHastingsRatio;
		
		for (int seqSite = 0; seqSite < pathTree.getSequences().get(0).getSequence().length; seqSite ++) {
			PupkoOneSite(pathTree, seqSite);
		}
		
		// NielsenOneSite
		// don't need to traverse tree (we can work on m_branches directly, since now we have internal states already)
		// TO-DO
		for (int branchNr = 0; branchNr < pathTree.getBranches().size(); branchNr++) {
			if((branchNr != rootNr) && (branchNr != sudoRootNr)) {
				NielsenSampleOneBranch(pathTree, branchNr);
			}
		}
		
		pathTree.showSequences();
	    System.out.println("--------------------------------------------------------------------------");

		/*
		pathTree.setDummySeqInternalNodesAll();
				*/
		
		// if this is the first run, we need to come up with oldPathDensity first
		
		
		//
		
		// change the tree and set the tree to be dirty
		pathTree.setSomethingIsDirty(true);
		
		// MCMC to test our implementation:
		// case 1: reject
		
		// case 2: accept
		
		// calculate new path density
		
		// 1. Pupko: propose ancestral states
		
		// 2. Nielsen: propose substitution events
		
		// combine Pupko and Nielsen's partial density
		
		// calculate hastings ratio
		newPathDensity = 0.1;
		oldPathDensity = 0.2;
		fHastingsRatio = newPathDensity / oldPathDensity;
		// to make sure MCMC will reject everytime
		fHastingsRatio = Double.NEGATIVE_INFINITY;
		
		// to make sure MCMC will accept everytime
		// fHastingsRatio = Double.POSITIVE_INFINITY;
		
		return fHastingsRatio;
	}
	
	public List<double[]> PupkoOneSite(PathTree pathTree, int seqSite){		
		int rootNr = pathTree.getRoot().getNr();
		int sudoRootNr = 0;
		for (Node childNode : pathTree.getRoot().getChildren()) {
			if(!childNode.isLeaf()) {
				sudoRootNr = childNode.getNr();
			}
		}
		
		// if the seq state at this site of the current node "N" is "j" and its parent seq state is "i"
		// then P(N=j | N_parent = i) is the "i*4 + j"th item in the "N"th Matrix of the list
		List<double[]> pMatrices = new ArrayList<double []>();
		for (int i=0; i< pathTree.getNodeCount(); i++){
			pMatrices.add(new double[4*4]);
		}

		
		
		/* deal with leaf nodes first (since they can be readily calculated) */
		for (int leafNodeNr = 0; leafNodeNr < pathTree.getLeafNodeCount(); leafNodeNr++) {
			setLeafpMatrix(pMatrices.get(leafNodeNr), pathTree.getNode(leafNodeNr).getLength(), pathTree.getSequences().get(leafNodeNr).getSequence()[seqSite]);
		}
		
		// create a list to store numbers for internal nodes
		List<Integer> internalNodeIndices = new ArrayList<Integer>();
		for (int internalNodeIndex = pathTree.getLeafNodeCount(); internalNodeIndex < pathTree.getNodeCount(); internalNodeIndex++) {
			internalNodeIndices.add(internalNodeIndex);
		}
		// remove the root number from internalNodeNum
		internalNodeIndices.remove(internalNodeIndices.indexOf(pathTree.getRoot().getNr()));

		while(internalNodeIndices.size() > 0 ){
			for (int i = 0; i < internalNodeIndices.size(); i++) {
				int currentNodeNr = internalNodeIndices.get(i);
				Node currentNode = pathTree.getNode(currentNodeNr);
				System.out.println("internalnode:" + currentNodeNr);
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
				
				if (currentNode.getLength() == 0){
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
						setSudoRoot(pathTree, currentNode.getNr(), seqSite, pMatrices.get(leftNum), pMatrices.get(rightNum), pMatrices.get(tipNodeNr));
						// trace back and set all internal nodes based on pMatrices
						setInternalNodes(pathTree, currentNode.getNr(), seqSite, pMatrices);
						internalNodeIndices.remove(i);
					}
				}
				
			}
		}

		// let the root have the same sequence as the sudoRoot
		MutableSequence newRootSeq = pathTree.getSequences().get(sudoRootNr).copy();
		pathTree.getSequences().get(rootNr).setSequence(newRootSeq.getSequence());
		
		return pMatrices;
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
	
	public double setSudoRoot(PathTree pathTree, int sudoRootNr, int seqSite, double [] leftpMatrix, double [] rightpMatrix, double[] tipNodepMatrix){
		
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
		
		System.out.println("CDF:" + Arrays.toString(nucleoStateCDF));
		
		// pick a nucleoState randomly based on nucleoStateProbs
		int sudoRootNucleoState = Randomizer.randomChoice(nucleoStateCDF);
		sudoRootNucleoStateProb = nucleoStateProbs[sudoRootNucleoState];
		pathTree.getSequences().get(sudoRootNr).getSequence()[seqSite] = sudoRootNucleoState;
		
		return sudoRootNucleoStateProb;
	}
	
	public void setInternalNodes(PathTree pathTree, int sudoRootNr, int seqSite, List<double []> pMatrices) {
		
		// get sudoRoot state
		int sudoRootNucleoState = pathTree.getSequences().get(sudoRootNr).getSequence()[seqSite];
		Node leftNode = pathTree.getNode(sudoRootNr).getLeft();
		Node rightNode = pathTree.getNode(sudoRootNr).getRight();
		traverseInternalNode(leftNode, pathTree, seqSite, pMatrices);
		traverseInternalNode(rightNode, pathTree, seqSite, pMatrices);
		// we need a recursive function here to traverse the tree and set states for internal nodes
		// idea:
		// traverse(left)
		// traverse(right)
		// traverse{
		//	traverse(left)
		//  traverse(right)
		// }
		//if() {
			
		//}
		
	}
	
	public void traverseInternalNode(Node node, PathTree pathTree, int seqSite, List<double []> pMatrices) {
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
			pathTree.getSequences().get(node.getNr()).getSequence()[seqSite] = thisNodeState;
			
			//recursive part:
			traverseInternalNode(node.getLeft() , pathTree, seqSite, pMatrices);
			traverseInternalNode(node.getRight() , pathTree, seqSite, pMatrices);
			
			// think about how to return value
			// return thisNodeStateProb;
		}
	}
	
	// matrix calculation helper function
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
	
	public double NielsenSampleOneBranch(PathTree pathTree, int branchNr) {
		return 0;
	}
}
