package evoprotein.proteinstructure;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import beast.core.Plugin;
import beast.util.Randomizer;

public class StructureEnv extends Plugin {
	
	List<double [][]> structEnv;
	
	// first row represents row marginal, second row represents column marginal
	List<double [][]> marginalProbMatrices;
	
	// pre-compute log would speed up my application
	List<double [][]> logStructEnv;
	List<double [][]> marginalLogProbMatrices;
	
	// initiate and validate
	public void initAndValidate(){
		parseStructureEnv();
	}
	
	public void parseStructureEnv(){
		// for now, generate matrices
		structEnv = mockUpMatrices(10);
		marginalProbMatrices = getMarginalProbMatrices();
		
		// tmp, used to test whether pre-compute log would speed up my application
		marginalLogProbMatrices = calculateMarginalLogProbMatrices(marginalProbMatrices);
		// after computing marginalLogProbMatrices, we can subtract marginalLogProbs away from logStructEnv
		logStructEnv = calculateLogStructEnv(structEnv, marginalLogProbMatrices);
		
	}
	
	List<double [][]>  calculateMarginalLogProbMatrices(List<double [][]> marginalProbMatrices){
		List<double [][]> marginalLogProbMatrices  = new ArrayList<double [][]>();
		
		for (int envNum = 0; envNum < marginalProbMatrices.size(); envNum++) {
			marginalLogProbMatrices.add(getLogMatrix(marginalProbMatrices.get(envNum)));
		}
		
		return marginalLogProbMatrices;
	}
	
	List<double [][]> calculateLogStructEnv(List<double [][]> structEnv, List<double [][]> marginalLogProbMatrices){
		List<double [][]> logStructEnv  = new ArrayList<double [][]>();
		// take each env and pre-compute log value for each matrix
		for (int envNum = 0; envNum < structEnv.size(); envNum++) {
			logStructEnv.add(getLogMatrixAdjustedByMarginal(structEnv.get(envNum), marginalLogProbMatrices.get(envNum)));
		}
		return logStructEnv;
	}
	
	double [][] getLogMatrix(double [][] matrix){
		int rowDim = matrix.length;
		int colDim = matrix[0].length;
		double [][] logMatrix = new double[rowDim][colDim];
		
		for(int i = 0; i < rowDim; i++){
			for (int j = 0 ;j< colDim; j++){
				logMatrix[i][j] = Math.log(matrix[i][j]);
			}
		}
		
		return logMatrix;
	}
	
	double [][] getLogMatrixAdjustedByMarginal(double [][] matrix, double[][] marginalLogMatrix){
		int rowDim = matrix.length;
		int colDim = matrix[0].length;
		double [][] logMatrix = new double[rowDim][colDim];
		
		for(int i = 0; i < rowDim; i++){
			for (int j = 0 ;j< colDim; j++){
				logMatrix[i][j] = Math.log(matrix[i][j]) - marginalLogMatrix[0][i] - marginalLogMatrix[1][j];
			}
		}
		
		return logMatrix;
	}	
	
	// getter
	public List<double[][]> getLogStructEnv(){
		return logStructEnv;
	}
	
	List<double [][]> getMarginalProbMatrices () {
		List<double [][]> marginalMatrices = new ArrayList<double [][]>();
		for (int envNum = 0; envNum < structEnv.size(); envNum++) {
			marginalMatrices.add(getMarginalProbMatrix(structEnv.get(envNum)));
		}
		return marginalMatrices;
	}
	
	public double [][] getMarginalProbMatrix(double[][] structEnvMatrix){
		int structEnvMatrixDim = structEnvMatrix.length;
		// create [2][61] matrix
		double [][] marginalMatrix = new double [2][structEnvMatrixDim];
		// calculate row & col marginal
		for (int rowNr = 0 ; rowNr < structEnvMatrixDim; rowNr++) {
			marginalMatrix[0][rowNr] = getRowSum(rowNr, structEnvMatrix);
			int colNr = rowNr;
			marginalMatrix[1][colNr] = getColSum(colNr, structEnvMatrix);
		}
		// return result
		return marginalMatrix;
	}
	
	double getRowSum (int rowNr, double[][] structEnvMatrix) {
		double rowSum = 0;
		for (int colNr = 0; colNr < structEnvMatrix[0].length; colNr++) {
			rowSum += structEnvMatrix[rowNr][colNr];
		}
		return rowSum;
	}
	
	double getColSum (int colNr, double[][] structEnvMatrix) {
		double colSum = 0;
		for (int rowNr = 0; rowNr < structEnvMatrix.length; rowNr++) {
			colSum += structEnvMatrix[rowNr][colNr];
		}
		return colSum;
	}
	
	public int getStructEnvNum(){
		return structEnv.size();
	}
	
	// needs efficiency boost (pre-calculation can do the trick)
	public double getProb(int structEnvNumber, int firstCodonType, int secondCodonType){
		double prob = structEnv.get(structEnvNumber)[firstCodonType][secondCodonType];
		double marginalProb = getFirstCodonMarginalProb(structEnvNumber, firstCodonType) * getSecondCodonMarginalProb(structEnvNumber, secondCodonType);
		prob = prob / marginalProb;
		return prob;
	}

	// already taken marginal prob into account
	public double getLogProb(int structEnvNumber, int firstCodonType, int secondCodonType){
		return logStructEnv.get(structEnvNumber)[firstCodonType][secondCodonType];
	}
	
	/*
	public double getLogProb(int structEnvNumber, int firstCodonType, int secondCodonType){
		double logProb = logStructEnv.get(structEnvNumber)[firstCodonType][secondCodonType];
		double marginalLogProb = getFirstCodonMarginalLogProb(structEnvNumber, firstCodonType) + getSecondCodonMarginalLogProb(structEnvNumber, secondCodonType);
		logProb = logProb - marginalLogProb;
		return logProb;
	}
	*/
	
	public double getFirstCodonMarginalProb(int structEnvNumber, int firstCodonType){
		return marginalProbMatrices.get(structEnvNumber)[0][firstCodonType];
	}
	
	public double getSecondCodonMarginalProb(int structEnvNumber, int secondCodonType){
		return marginalProbMatrices.get(structEnvNumber)[1][secondCodonType];
	}
	
	public double getFirstCodonMarginalLogProb(int structEnvNumber, int firstCodonType){
		return marginalLogProbMatrices.get(structEnvNumber)[0][firstCodonType];
	}
	
	public double getSecondCodonMarginalLogProb(int structEnvNumber, int secondCodonType){
		return marginalLogProbMatrices.get(structEnvNumber)[1][secondCodonType];
	}
	
	/*
	public double getFirstCodonMarginalProb(int structEnvNumber, int firstCodonType){
		double firstCodonMarginalProb = 0;
		int rowNr = firstCodonType;
		for (int i = 0 ; i < structEnv.get(0)[0].length ; i++) {
			firstCodonMarginalProb += structEnv.get(structEnvNumber)[rowNr][i];
		}
		return firstCodonMarginalProb;
	}
	
	public double getSecondCodonMarginalProb(int structEnvNumber, int secondCodonType){
		double secondCodonMarginalProb = 0;
		int colNr = secondCodonType;
		for (int i = 0; i < structEnv.get(0).length ; i++) {
			secondCodonMarginalProb += structEnv.get(structEnvNumber)[i][colNr];
		}
		return secondCodonMarginalProb;
	}
	*/
	
	// for testing purpose only
	public List<double [][]> mockUpMatrices(int numberOfMatrices){
		List<double [][]> matrices = new ArrayList<double [][]>();
		for (int i = 0 ; i < numberOfMatrices; i++) {
			matrices.add(mockUpMatrix());
		}
		return matrices;
	}
	
	public double[][] mockUpMatrix(){

		double [] codonFreq = new double[]{0.00983798,0.01745548,0.00222048,0.01443315,
				0.00844604,0.01498576,0.00190632,0.01239105,
				0.01064012,0.01887870,
				0.00469486,0.00833007,0.00688776,
				0.01592816,0.02826125,0.00359507,0.02336796,
				0.01367453,0.02426265,0.00308642,0.02006170,
				0.01722686,0.03056552,0.00388819,0.02527326,
				0.00760121,0.01348678,0.00171563,0.01115161,
				0.01574077,0.02792876,0.00355278,0.02309304,
				0.01351366,0.02397721,0.00305010,0.01982568,
				0.01702419,0.03020593,0.00384245,0.02497593,
				0.00751178,0.01332811,0.00169545,0.01102042,
				0.02525082,0.04480239,0.00569924,0.03704508,
				0.02167816,0.03846344,0.00489288,0.03180369,
				0.02730964,0.04845534,0.00616393,0.04006555,
				0.01205015,0.02138052,0.00271978,0.01767859};
		
		double [][] mockupmatrix = new double [61][61];

		for (int i = 0; i < 61; i++) {
			for (int j = 0; j < 61; j++) {
				mockupmatrix[i][j] = codonFreq[i] * codonFreq[j];
			}
		}
		
		return mockupmatrix;
	}
	
}
