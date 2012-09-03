package beast.evolution.tree;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.datatype.Nucleotide;


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
	
	// setters
	
	// for testing purpose only
	public void setDummySeqInternalNodes (int [] dummySeq){
		for(int i=0; i<nodeCount; i++){
			m_sequences[i] = dummySeq;
		}
	}
	
	public void showSequences() {
		Nucleotide nucleo = new Nucleotide();
		for (int i = 0; i < m_nodes.length; i++) {
			int[] tmp_seq = m_sequences[i];
			/**
			 * System.out.println("The length of the sequence is:" +
			 * aligned_seq.getSiteCount());
			 */
			String real_seq = nucleo.state2string(tmp_seq);
			System.out.println("The node is:" + i);
			if (m_nodes[i].isLeaf()) {
				System.out.println("The ID is:" + m_nodes[i].getID()
						+ " and length is: " + m_nodes[i].getLength());
			}
			System.out.println(" and the Seq is: \n" + real_seq);
		}
	}
	
	// getters
	public int[][] getSequences(){
		return m_sequences;
	}
}
