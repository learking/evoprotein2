package beast.evolution.substitutionmodel;

import beast.core.CalculationNode;
import beast.evolution.datatype.DataType;
import beast.evolution.tree.Node;

public class ProteinCodingDNASubstModel extends CalculationNode implements SubstitutionModel {

	// determine inputs
	
	// define variables here
	
	// initAndValidate
	
	
	
	@Override
	public void getTransitionProbabilities(Node node, double fStartTime,
			double fEndTime, double fRate, double[] matrix) {
		// TODO Auto-generated method stub
		
	}

	// most important implementation
	@Override
	public double[] getRateMatrix(Node node) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double[] getFrequencies() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getStateCount() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public EigenDecomposition getEigenDecomposition(Node node) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean canReturnComplexDiagonalization() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean canHandleDataType(DataType dataType) throws Exception {
		// TODO Auto-generated method stub
		return false;
	}
	
}
