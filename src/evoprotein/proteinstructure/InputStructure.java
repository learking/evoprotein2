package evoprotein.proteinstructure;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.util.Randomizer;

public class InputStructure extends BEASTObject {
	
	// for debugging only
	int simulatedCodonNr = 137;
	
	public Input<StructureEnv> structureEnv = new Input<StructureEnv>("structureEnv", "Structure Environment", Validate.REQUIRED);
	public Input<SolventAccessibility> solventAccessibility = new Input<SolventAccessibility>("solventAccessibility", "solvent accessibility categories", Validate.REQUIRED);
	
	public Input<String> m_fistOrderTermsFile = new Input<String>("firstOrderTerms", "firstOrderTerm file", Validate.REQUIRED);// not required, Beauti may need this for example    
	public Input<String> m_interactionTermsFile = new Input<String>("interactionTerms", "interactionTerm file", Validate.REQUIRED);// not required, Beauti may need this for example
	
	
	protected int [] firstOrderTerms; 
	
	// need a sparse matrix here to store pairwise info, however, since the # of entries is much smaller than ...
	protected int [][] interactionTerm2EnvMap;
	
	
	// introduce this in order to run faster, but does it help?
	double[][] solventCategories;
	List<double[][]> logStructEnv;
	
	// initiate and validate
	public void initAndValidate() throws IOException{
		System.out.println(m_interactionTermsFile.get());
		parseInputStructure();
		
		// try to speed up
		solventCategories = solventAccessibility.get().getSolventCategories(); 
		logStructEnv = structureEnv.get().getLogStructEnv();
	}
	
	// file parser
	public void parseInputStructure() throws IOException{
		
		// set p(i|c_sol)
		parseFirstOrderTerms();
		
		// set p(m, n | c_sol, c_distance)
		parseInteractionTerm2EnvMap();
	}
	
	// first order
	public void parseFirstOrderTerms() throws IOException{
		// for now, mock up a list that logs which solvent category this position belongs to
		// length of the list should be equal to AAseqLength
		//firstOrderTerms = mockUpFirstOrderTerms();
		firstOrderTerms = readFirstOrderTerms();
	}
	
	// second order
	public void parseInteractionTerm2EnvMap() throws IOException{
		// sparse matrix might be a better choice, but for now, I will use matrix
		// for now, mock up a list that logs which environment this pair of AA under consideration belongs to
		//interactionTerm2EnvMap = mockUpInteractionTerm2EnvMap();
		interactionTerm2EnvMap = readInteractionTerm2EnvMap();
	}
	
	// getters
	
	public double getFirstOrderLogProb(int codonPosition, int codonType) {
		int firstOrderTermCategory = firstOrderTerms[codonPosition];
		return solventCategories[firstOrderTermCategory][codonType];
	}
	
	// needs efficiency boost
	public double getInteractionLogProb(int firstCodonPosition, int secondCodonPosition, int firstCodonType, int secondCodonType){
		int structEnvNumber = interactionTerm2EnvMap[firstCodonPosition][secondCodonPosition];
		return logStructEnv.get(structEnvNumber)[firstCodonType][secondCodonType];
	}
	
	
	
	// read in terms
	int[] readFirstOrderTerms() throws IOException {
		FileReader reader = new FileReader(m_fistOrderTermsFile.get());
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
	
	// read interactionTerms
	int[][] readInteractionTerm2EnvMap() throws IOException{
		FileReader reader = new FileReader(m_interactionTermsFile.get());
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
	
	// remove deletion-caused-gap related first order and second order terms
	public void removeGapRelatedTerms(Set<Integer> deletionPositions) throws Exception{
		int newDim = getFirstOrderDim() - deletionPositions.size(); 
		int[] newFirstOrderTerms = new int[newDim];
		int[][] newInteractionTerm2EnvMap = new int[newDim][newDim];
		
		// delete first order terms by creating a new int[] and replace the original one with the new one
		int j  = 0;
		for(int i = 0; i < getFirstOrderDim(); i++){
			// add to new int[] only when it is not deletion site
			if(!deletionPositions.contains(i*3)){
				newFirstOrderTerms[j] = getFirstOrderTerm(i);
				
				j++;
				if(j > newDim){
					throw new Exception("J shouldn't exceed newDim");
				}
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
						newInteractionTerm2EnvMap[m][n] = getInteractionTerm(rowNr, colNr);
					
						// update after operation
						n++;
					}
				}
						
				m++;		
			}
		
		}
		
		updateFirstOrderTerms(newFirstOrderTerms);
		updateInteractionTerms(newInteractionTerm2EnvMap);
		
	}
	
	void updateFirstOrderTerms(int[] newFirstOrderTerms){
		firstOrderTerms = newFirstOrderTerms;
	}
	
	void updateInteractionTerms(int[][] newInteractionTerm2EnvMap){
		interactionTerm2EnvMap = newInteractionTerm2EnvMap;
	}
	
	// getter
	int getFirstOrderDim(){
		return firstOrderTerms.length;
	}
	
	int getInteractionDim(){
		return interactionTerm2EnvMap.length;
	}
	
	int getFirstOrderTerm(int i){
		return firstOrderTerms[i];
	}
	
	int getInteractionTerm(int rowNr, int colNr){
		return interactionTerm2EnvMap[rowNr][colNr];
	}
	
	// for testing only
	// for test only
	
	public void writeFirstOrder() throws IOException {
		// mac
		//FileWriter fstream = new FileWriter("/Users/kwang2/Desktop/firstOrderTerms.txt");
		// linux
		FileWriter fstream = new FileWriter("/home/kuangyu/Desktop/firstOrderTerms.txt");
		BufferedWriter out = new BufferedWriter(fstream);
		  
			int [] firstOrder = new int [simulatedCodonNr];
			double [] pdf = new double [10];
			
			double tempSum = 0;
			for (int i = 0; i < pdf.length; i++) {
				tempSum += 1.0 / (double) pdf.length;
				pdf[i] = tempSum;
			}
			
			for (int i = 0; i < firstOrder.length; i++) {
				firstOrder[i] = Randomizer.randomChoice(pdf);
			}
			
			// print the entire vector
			for(int i=0; i < firstOrder.length; i++) {
				out.write(firstOrder[i] + " ");
			}
			out.write("\n");
			out.write("//");
		  
		  //Close the output stream
		  out.close();
	}
	
	public void printFirstOrder(){
		// for categories available, randomly assign category to each position
		int [] firstOrder = new int [simulatedCodonNr];
		double [] pdf = new double [10];
		
		double tempSum = 0;
		for (int i = 0; i < pdf.length; i++) {
			tempSum += 1.0 / (double) pdf.length;
			pdf[i] = tempSum;
		}
		
		for (int i = 0; i < firstOrder.length; i++) {
			firstOrder[i] = Randomizer.randomChoice(pdf);
		}
		
		// print the entire vector
		for(int i=0; i < firstOrder.length; i++) {
			System.out.print(firstOrder[i] + " ");
		}
		System.out.println();
	}
	
	public void writeInteraction2EnvMap() throws IOException {
		//mac
		  //FileWriter fstream = new FileWriter("/Users/kwang2/Desktop/interactionTerms.txt");
		//linux
		FileWriter fstream = new FileWriter("/home/kuangyu/Desktop/interactionTerms.txt");
		BufferedWriter out = new BufferedWriter(fstream);
		  
			int[][] interactionTerm2EnvMap = new int[simulatedCodonNr][simulatedCodonNr];

			int structEnvNum = 10;
			
			double tempSum = 0;
			double [] pdf = new double [structEnvNum];
			for (int i = 0; i < pdf.length; i++) {
				tempSum += 1.0 / (double) pdf.length;
				pdf[i] = tempSum;
			}
			
			for (int rowNr = 0 ; rowNr < interactionTerm2EnvMap.length; rowNr++) {
				for (int colNr = rowNr; colNr < interactionTerm2EnvMap[0].length; colNr++) {
					interactionTerm2EnvMap[rowNr][colNr] = Randomizer.randomChoice(pdf);
				}
			}
			
			// print the entire matrix
			for (int rowNr = 0 ; rowNr < interactionTerm2EnvMap.length; rowNr++) {
				for (int colNr = 0; colNr < interactionTerm2EnvMap[0].length; colNr++) {
					out.write(interactionTerm2EnvMap[rowNr][colNr] + " ");
				}
				out.write("\n");
			}
		  out.write("//");
		  //Close the output stream
		  out.close();
	}
	
	public void printInteraction2EnvMap(){
		int[][] interactionTerm2EnvMap = new int[simulatedCodonNr][simulatedCodonNr];

		int structEnvNum = 10;
		
		double tempSum = 0;
		double [] pdf = new double [structEnvNum];
		for (int i = 0; i < pdf.length; i++) {
			tempSum += 1.0 / (double) pdf.length;
			pdf[i] = tempSum;
		}
		
		for (int rowNr = 0 ; rowNr < interactionTerm2EnvMap.length; rowNr++) {
			for (int colNr = rowNr; colNr < interactionTerm2EnvMap[0].length; colNr++) {
				interactionTerm2EnvMap[rowNr][colNr] = Randomizer.randomChoice(pdf);
			}
		}
		
		// print the entire matrix
		for (int rowNr = 0 ; rowNr < interactionTerm2EnvMap.length; rowNr++) {
			for (int colNr = 0; colNr < interactionTerm2EnvMap[0].length; colNr++) {
				System.out.print(interactionTerm2EnvMap[rowNr][colNr] + " ");
			}
			System.out.println();
		}
		
	
	}
	
	public void showFirstOrderTerms() {
		// print the entire vector
		for(int i=0; i < firstOrderTerms.length; i++) {
			System.out.print(firstOrderTerms[i] + " ");
		}
		System.out.println();
	}
	
	public void showInteractions() {
		for (int rowNr = 0 ; rowNr < interactionTerm2EnvMap.length; rowNr++) {
			for (int colNr = 0; colNr < interactionTerm2EnvMap[0].length; colNr++) {
				System.out.print(interactionTerm2EnvMap[rowNr][colNr] + " ");
			}
			System.out.println();
		}
	}
}
