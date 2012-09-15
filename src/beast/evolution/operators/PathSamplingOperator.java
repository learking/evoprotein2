/**
 * 
 */
package beast.evolution.operators;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.PathTree;

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

		// Is there any way to utilize old path density without recalculation?
		// for example, if the proposal gets rejected
		double oldPathDensity, newPathDensity, fHastingsRatio;
		
		List<double []> pMatrices = PupkoOneSite(pathTree, 0);
		
		for (double[] tmpMatrix : pMatrices){
			System.out.println(Arrays.toString(tmpMatrix));
		}
		
		/*
		pathTree.setDummySeqInternalNodesAll();
		pathTree.showSequences();
	    System.out.println("--------------------------------------------------------------------------");
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
		final double [] zeroArray = new double [4*4];
		
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
					if ((!pMatrices.get(leftNum).equals(zeroArray))&&(!pMatrices.get(rightNum).equals(zeroArray))) {
						// set up pMatrix for this internal node
						setInternalpMatrix(pMatrices.get(currentNodeNr), pMatrices.get(leftNum), pMatrices.get(rightNum), currentNode.getLength());
						internalNodeIndices.remove(i);
					}
				}
				
				if (currentNode.getLength() == 0){
					int leftNum = currentNode.getLeft().getNr();
					int rightNum = currentNode.getRight().getNr();
				    // if pMatrices of both left and right child node has been calculated
					if ((!pMatrices.get(leftNum).equals(zeroArray))&&(!pMatrices.get(rightNum).equals(zeroArray))) {
						// set up pMatrix for this internal node
						//setSudoRootpMatrix(pMatrices.get(currentNodeNr), pMatrices.get(leftNum), pMatrices.get(rightNum), currentNode.getLength());
						internalNodeIndices.remove(i);
					}
				}
				
			}
		}

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
	
	public double setSudoRoot(double [] leftpMatrix, double [] rightpMatrix, double[] tipNodepMatrix){
		double sudoRootNucleoStateProb = 0.0;
		final double [] frequences = m_substitutionModel.getFrequencies();
		double [] nucleoStateProbs = new double [4*4];
		for(int nucleoState=0; nucleoState < 4; nucleoState++){
			double nucleoStateFreq = frequences[nucleoState];
			double leftpMatrixRowSum = getpMatrixRowSum(leftpMatrix, nucleoState);
			double rightpMatrixRowSum = getpMatrixRowSum(rightpMatrix, nucleoState);
			double tipNodepMatrixRowSum = getpMatrixRowSum(tipNodepMatrix, nucleoState);
			nucleoStateProbs[nucleoState] = nucleoStateFreq * leftpMatrixRowSum * rightpMatrixRowSum * tipNodepMatrixRowSum;
		}
		// TO-DO: select a state based on nucleoStateProbs
		
		return sudoRootNucleoStateProb;
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
	
}
