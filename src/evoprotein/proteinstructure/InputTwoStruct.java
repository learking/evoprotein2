package evoprotein.proteinstructure;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

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
	
	// initiate and validate
	public void initAndValidate() throws IOException{
		parseInputTwoStruct();
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
		return solventAccessibility.get().getLogProb(firstOrderTermCategory, codonType); 
	}
	
	public double getInteractionLogProb(int firstCodonPosition, int secondCodonPosition, int firstCodonType, int secondCodonType, int[][] interactionTerm2EnvMap) {
		int structEnvNumber = interactionTerm2EnvMap[firstCodonPosition][secondCodonPosition];
		return structureEnv.get().getLogProb(structEnvNumber, firstCodonType, secondCodonType);
	}
	
	/*
	public double getInteractionLogProb(int firstCodonPosition, int secondCodonPosition, int firstCodonType, int secondCodonType){
		int structEnvNumber = interactionTerm2EnvMap[firstCodonPosition][secondCodonPosition];
		return structureEnv.get().getLogProb(structEnvNumber, firstCodonType, secondCodonType);
	}
	*/
	
	public double getFirstOrderRatio(int codonDifferPosition, int differCodonType, int originalCodonType, double fNow){
		double firstOrderRatio = 0;
		double firstOrderStructARatio = 0;
		double firstOrderStructBRatio = 0;
		
		firstOrderStructARatio = getFirstOrderLogProb(codonDifferPosition, differCodonType, firstOrderStructA) - getFirstOrderLogProb(codonDifferPosition, originalCodonType, firstOrderStructA);
		firstOrderStructBRatio = getFirstOrderLogProb(codonDifferPosition, differCodonType, firstOrderStructB) - getFirstOrderLogProb(codonDifferPosition, originalCodonType, firstOrderStructB);
	    firstOrderRatio = fNow*firstOrderStructARatio + (1-fNow)*firstOrderStructBRatio;
		return firstOrderRatio;
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
}
