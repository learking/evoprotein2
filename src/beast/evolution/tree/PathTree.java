package beast.evolution.tree;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;


@Description("Tree that enables path sampling")
public class PathTree extends Tree {
	
	public Input<Alignment> m_alignment = new Input<Alignment>("alignment",
			"alignment that contains sequence data");
	
	int[][] m_sequences;
	
	@Override
	public void initAndValidate() throws Exception {
			super.initAndValidate();
			m_sequences = new int[nodeCount][m_alignment.get().getSiteCount()];
	}
	
	
	
	//getters
	public int[][] getSequences(){
		return m_sequences;
	}
}
