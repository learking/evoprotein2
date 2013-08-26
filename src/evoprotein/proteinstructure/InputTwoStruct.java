package evoprotein.proteinstructure;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import beast.core.Input;
import beast.core.Plugin;
import beast.core.Input.Validate;

public class InputTwoStruct extends Plugin {

	public Input<StructureEnv> structureEnv = new Input<StructureEnv>("structureEnv", "Structure Environment", Validate.REQUIRED);
	public Input<SolventAccessibility> solventAccessibility = new Input<SolventAccessibility>("solventAccessibility", "solvent accessibility categories", Validate.REQUIRED);

	// unlike previously, now we have four input files specifying both first order and interaction terms for struct A and B
	public Input<String> m_fistOrderStructA= new Input<String>("firstOrderStructA", "firstOrderTerm file for struct A", Validate.REQUIRED);// not required, Beauti may need this for example    
	public Input<String> m_interactionTermsStructA = new Input<String>("interactionTermsStructA", "interactionTerm file for struct A", Validate.REQUIRED);// not required, Beauti may need this for example
	public Input<String> m_fistOrderStructB= new Input<String>("firstOrderStructB", "firstOrderTerm file for struct B", Validate.REQUIRED);// not required, Beauti may need this for example    
	public Input<String> m_interactionTermsStructB = new Input<String>("interactionTermsStructB", "interactionTerm file for struct B", Validate.REQUIRED);// not required, Beauti may need this for example
		
	int [] firstOrderStructA;
	int [] firstOrderStructB;	
	int [][] interactionTermsStructA2EnvMap;
	int [][] interactionTermsStructB2EnvMap;
	
	// introduce this in order to run faster, but does it help?
	double[][] solventCategories;
	List<double[][]> logStructEnv;
	
	// initiate and validate
	public void initAndValidate() throws IOException{
		parseInputTwoStruct();
		
		// try to speed up
		solventCategories = solventAccessibility.get().getSolventCategories(); 
		logStructEnv = structureEnv.get().getLogStructEnv();
	}
	
	public void parseInputTwoStruct() throws IOException{
		// set p(i|c_sol)
		parseFirstOrderTerms();
		// set p(m, n | c_sol, c_distance)
		parseInteractionTerm2EnvMap();
	}
	
	void parseFirstOrderTerms() throws IOException{
		// parse first order for struct A
		firstOrderStructA = readFirstOrderTerms(m_fistOrderStructA.get());
		// parse first order for struct B
		firstOrderStructB = readFirstOrderTerms(m_fistOrderStructB.get());
	}
	
	// read in terms
	int[] readFirstOrderTerms(String firstOrderFile) throws IOException {
		FileReader reader = new FileReader(firstOrderFile);
		BufferedReader bufferedReader = new BufferedReader(reader);

		List<String> lines = new ArrayList<String>();
		String line = "";
		while(!(line = bufferedReader.readLine()).trim().equals("//")) {
			lines.add(line);
		}
		
		int vectorDim = lines.get(0).split("\\s+").length;
		int[] vector = new int[vectorDim];
		
		String[] firstOrderTermsStr = lines.get(0).split("\\s+");
		for (int item = 0; item < vectorDim; item++) {
			vector[item] = Integer.parseInt(firstOrderTermsStr[item]);
		}
		return vector;
	}
	
	// 
	void parseInteractionTerm2EnvMap() throws IOException{
		// parse interactions for struct A
		interactionTermsStructA2EnvMap = readInteractionTerm2EnvMap(m_interactionTermsStructA.get());
		// parse interactions for struct B	
		interactionTermsStructB2EnvMap = readInteractionTerm2EnvMap(m_interactionTermsStructB.get());
	}
	
	int[][] readInteractionTerm2EnvMap(String interactionTerms2EnvMapFile) throws IOException{
		FileReader reader = new FileReader(interactionTerms2EnvMapFile);
		BufferedReader bufferedReader = new BufferedReader(reader);
		
		List<String> lines = new ArrayList<String>();
		String line = "";
		while(!(line = bufferedReader.readLine()).trim().equals("//")) {
			lines.add(line);
		}
		
		int matrixDim = lines.get(0).split("\\s+").length;
		int[][] matrix = new int[matrixDim][matrixDim];
		for (int i = 0; i < lines.size(); i++) {
			String[] matrixEntriesStr = lines.get(i).split("\\s+");
			int[] matrixEntries = new int[matrixEntriesStr.length];
			for (int item = 0; item < matrixEntries.length; item++) {
				matrixEntries[item] = Integer.parseInt(matrixEntriesStr[item]);
				//System.out.print(matrixEntries[item] + "|");
			}
			matrix[i] = matrixEntries;
		}
		
		return matrix;
	}
	
	// getter
	
	public double getFirstOrderLogProb(int codonPosition, int codonType, int[] firstOrderTerms) {
		int firstOrderTermCategory = firstOrderTerms[codonPosition];
		
		// after introducing solventAccessibility		
		return solventCategories[firstOrderTermCategory][codonType];
	}
	
	public double getInteractionLogProb(int firstCodonPosition, int secondCodonPosition, int firstCodonType, int secondCodonType, int[][] interactionTerm2EnvMap) {
		int structEnvNumber = interactionTerm2EnvMap[firstCodonPosition][secondCodonPosition];

		// after introducing logStructEnv
		return logStructEnv.get(structEnvNumber)[firstCodonType][secondCodonType];
	}
	
	public double getFirstOrderRatio(int codonDifferPosition, int differCodonType, int originalCodonType, double fNow){
		double firstOrderRatio = 0;
		double firstOrderStructARatio = 0;
		double firstOrderStructBRatio = 0;
		
		firstOrderStructARatio = getFirstOrderLogProb(codonDifferPosition, differCodonType, firstOrderStructA) - getFirstOrderLogProb(codonDifferPosition, originalCodonType, firstOrderStructA);
		firstOrderStructBRatio = getFirstOrderLogProb(codonDifferPosition, differCodonType, firstOrderStructB) - getFirstOrderLogProb(codonDifferPosition, originalCodonType, firstOrderStructB);
	    firstOrderRatio = fNow*firstOrderStructARatio + (1-fNow)*firstOrderStructBRatio;
		return firstOrderRatio;
	}

	// instead of calculating ratio, this function calculates firstOrderLogP for a given sequence at a given position
	public double getFirstOrderSeqLogP(int codonType, int codonPosition, double fNow){
		
		double firstOrderSeqLogProbStructA = getFirstOrderLogProb(codonPosition, codonType, firstOrderStructA);
		double firstOrderSeqLogProbStructB = getFirstOrderLogProb(codonPosition, codonType, firstOrderStructB);
		
		return fNow*firstOrderSeqLogProbStructA + (1-fNow)*firstOrderSeqLogProbStructB;
	}
	
	public double getInteractionRatio(int[] codonArrayI, int leftBound, int rightBound, int codonDifferPosition, int differCodon, double fNow) {
		double interactionRatio = 0;
		double interactionStructARatio = 0;
		double interactionStructBRatio = 0;

		// calculate interaction ratio for both structs independently (or jointly?)
		for (int m = leftBound; m < codonDifferPosition; m++) {
			//interactionRatio += getInteractionLogProb(m, codonDifferPosition, codonArrayI[m], differCodon) - getInteractionLogProb(m, codonDifferPosition, codonArrayI[m], codonArrayI[codonDifferPosition]);
			interactionStructARatio = getInteractionLogProb(m, codonDifferPosition, codonArrayI[m], differCodon, interactionTermsStructA2EnvMap) 
					- getInteractionLogProb(m, codonDifferPosition, codonArrayI[m], codonArrayI[codonDifferPosition], interactionTermsStructA2EnvMap);
			
			interactionStructBRatio = getInteractionLogProb(m, codonDifferPosition, codonArrayI[m], differCodon, interactionTermsStructB2EnvMap) 
					- getInteractionLogProb(m, codonDifferPosition, codonArrayI[m], codonArrayI[codonDifferPosition], interactionTermsStructB2EnvMap);
		}	
		for (int n = codonDifferPosition + 1 ; n < rightBound; n++) {
			interactionStructARatio += getInteractionLogProb(codonDifferPosition, n, differCodon, codonArrayI[n], interactionTermsStructA2EnvMap) 
			- getInteractionLogProb(codonDifferPosition, n, codonArrayI[codonDifferPosition], codonArrayI[n], interactionTermsStructA2EnvMap);
			
			interactionStructBRatio += getInteractionLogProb(codonDifferPosition, n, differCodon, codonArrayI[n], interactionTermsStructB2EnvMap) 
			- getInteractionLogProb(codonDifferPosition, n, codonArrayI[codonDifferPosition], codonArrayI[n], interactionTermsStructB2EnvMap);
		}	
		
		interactionRatio = fNow*interactionStructARatio + (1-fNow)*interactionStructBRatio;
		return interactionRatio;
	}
	
	public double getRootSeqLogP(int[] rootCodonSeq, double fNow){
		double rootSeqLogP = 0;
		
    	for(int codonPosition = 0; codonPosition < rootCodonSeq.length; codonPosition++) {
    		double firstOrderLogProbA = getFirstOrderLogProb(codonPosition, rootCodonSeq[codonPosition], firstOrderStructA);
    		double firstOrderLogProbB = getFirstOrderLogProb(codonPosition, rootCodonSeq[codonPosition], firstOrderStructB);
    		rootSeqLogP += fNow*firstOrderLogProbA + (1-fNow)*firstOrderLogProbB;
    	}
    	
    	// deal with second order terms
    	for(int firstCodonPosition = 0; firstCodonPosition < (rootCodonSeq.length-1); firstCodonPosition++) {
    		for(int secondCodonPosition = (firstCodonPosition+1); secondCodonPosition < rootCodonSeq.length; secondCodonPosition++) {
    			double interactionLogProb_A = getInteractionLogProb(firstCodonPosition, secondCodonPosition, rootCodonSeq[firstCodonPosition], rootCodonSeq[secondCodonPosition], interactionTermsStructA2EnvMap);
    			double interactionLogProb_B = getInteractionLogProb(firstCodonPosition, secondCodonPosition, rootCodonSeq[firstCodonPosition], rootCodonSeq[secondCodonPosition], interactionTermsStructB2EnvMap);
    			rootSeqLogP += fNow*interactionLogProb_A + (1-fNow)*interactionLogProb_B;
    		
    		}   		
    	}
		
		return rootSeqLogP;
	}

	public void removeGapRelatedTerms(Set<Integer> deletionPositions) {
		int newDim = getFirstOrderDim() - deletionPositions.size(); 
		
		int[] newFirstOrderStructA = new int[newDim];
		int[] newFirstOrderStructB = new int[newDim];		
		int[][] newInteractionTermsStructA2EnvMap = new int[newDim][newDim];
		int[][] newInteractionTermsStructB2EnvMap = new int[newDim][newDim];
		
		int j  = 0;
		for(int i = 0; i < getFirstOrderDim(); i++){
			// add to new int[] only when it is not deletion site
			if(!deletionPositions.contains(i*3)){
				newFirstOrderStructA[j] = getFirstOrderTermStructA(i);
				newFirstOrderStructB[j] = getFirstOrderTermStructB(i);
				
				j++;
				/*
				if(j > newDim){
					throw new Exception("J shouldn't exceed newDim");
				}
				*/
			}
		}
		
		int m = 0;
		// delete interaction terms by creating a new int[][] and replace the original one with the new one
		for(int rowNr = 0; rowNr < getInteractionDim(); rowNr++){
		
			System.out.println("interaction dim:" + getInteractionDim());			
			if(!deletionPositions.contains(rowNr*3)){
				// init n to be zero
				int n = 0;
				for(int colNr = 0; colNr < getInteractionDim(); colNr++){		
					if(!deletionPositions.contains(colNr*3)){
						newInteractionTermsStructA2EnvMap[m][n] = getInteractionTermStructA(rowNr, colNr);
						newInteractionTermsStructB2EnvMap[m][n] = getInteractionTermStructB(rowNr, colNr);
						// update after operation
						n++;
					}
				}
						
				m++;	
			}
		
		}
		
		updateFirstOrderTerms(newFirstOrderStructA, newFirstOrderStructB);
		updateInteractionTerms(newInteractionTermsStructA2EnvMap, newInteractionTermsStructB2EnvMap);
		
	}
		
	void updateFirstOrderTerms(int[] newFirstOrderStructA, int[] newFirstOrderStructB){
		firstOrderStructA = newFirstOrderStructA;
		firstOrderStructB = newFirstOrderStructB;
	}
	
	void updateInteractionTerms(int[][] newInteractionTermsStructA2EnvMap, int[][] newInteractionTermsStructB2EnvMap){
		interactionTermsStructA2EnvMap = newInteractionTermsStructA2EnvMap;
		interactionTermsStructB2EnvMap = newInteractionTermsStructB2EnvMap;
	}
	
	int getFirstOrderDim(){
		return firstOrderStructA.length;
	}
	
	int getInteractionDim(){
		return interactionTermsStructA2EnvMap.length;
	}
	
	int getFirstOrderTermStructA(int i){
		return firstOrderStructA[i];
	}

	int getFirstOrderTermStructB(int i){
		return firstOrderStructB[i];
	}
	
	int getInteractionTermStructA(int rowNr, int colNr){
		return interactionTermsStructA2EnvMap[rowNr][colNr];
	}

	int getInteractionTermStructB(int rowNr, int colNr){
		return interactionTermsStructB2EnvMap[rowNr][colNr];
	}
	
}
