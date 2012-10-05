package evoprotein.proteinstructure;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import beast.core.Plugin;
import beast.util.Randomizer;

public class StructureEnv extends Plugin {
	
	List<double [][]> structEnv;
	
	// initiate and validate
	public void initAndValidate(){
		parseStructureEnv();
	}
	
	public void parseStructureEnv(){
		// for now, generate matrices
		structEnv = mockUpMatrices(10);
	}
	
	// for testing purpose only
	public List<double [][]> mockUpMatrices(int numberOfMatrices){
		List<double [][]> matrices = new ArrayList<double [][]>();
		for (int i = 0 ; i < numberOfMatrices; i++) {
			matrices.add(mockUpMatrix());
		}
		return matrices;
	}
	
	public double[][] mockUpMatrix(){
		double [][] mockupmatrix = new double [61][61];
		double [] randomNumbers = new double [61*61];
		for (int i = 0; i < 61*61; i++) {
			randomNumbers[i] = Randomizer.nextDouble();
		}
		double [] normalizedNumbers = Randomizer.getNormalized(randomNumbers);
		for (int i = 0; i < 61; i++) {
			for (int j = 0; j < 61; j++) {
				mockupmatrix[i][j] = normalizedNumbers[i*61 + j];
			}
		}
		return mockupmatrix;
	}
	
}
