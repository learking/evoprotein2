package beast.evolution.tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import evoprotein.evolution.datatype.MutableSequence;
import evoprotein.evolution.substitution.SubstitutionEvent;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.datatype.Nucleotide;


@Description("Tree that enables path sampling")
public class PathTree extends Tree {
	
	public Input<Alignment> m_alignment = new Input<Alignment>("alignment",
			"alignment that contains sequence data");
	
	List<MutableSequence> m_sequences = new ArrayList<MutableSequence>();
	List<MutableSequence> m_storedsequences = new ArrayList<MutableSequence>();
	
	List<PathBranch> m_branches = new ArrayList<PathBranch>();
	List<PathBranch> m_storedbranches = new ArrayList<PathBranch>();
	
	@Override
	public void initAndValidate() throws Exception {
			super.initAndValidate();
			
			int sequenceLength = m_alignment.get().getSiteCount();
			for(int i=0; i<nodeCount; i++) {
				m_sequences.add(new MutableSequence(sequenceLength));
				if(m_nodes[i].isRoot()){
					m_branches.add(new PathBranch());
				}
				else{
					m_branches.add(new PathBranch(sequenceLength, m_nodes[i].getParent().getNr(), i));
				}
			}
			
			initLeafSequences();
	}
	
	private void initLeafSequences() throws Exception{
		for(int i=0; i< leafNodeCount;i++){
			String sequenceID = getNode(i).getID();
			m_sequences.get(i).setSequence(getSequenceByID(sequenceID));
		}
	}
	
	public int[] getSequenceByID(String sequenceID) throws Exception {
		Nucleotide nucleo = new Nucleotide();
		List<Integer> sequenceTarget = new ArrayList<Integer>();
		for (int i = 0; i < m_alignment.get().getNrTaxa(); i++) {
			if (m_alignment.get().m_pSequences.get().get(i).m_sTaxon.get()
					.toString() == sequenceID) {
				sequenceTarget = m_alignment.get().m_pSequences.get().get(i)
						.getSequence(nucleo);
			}
		}
		int[] sequence = new int[sequenceTarget.size()];
		for (int i = 0; i < sequence.length; i++) {
			sequence[i] = sequenceTarget.get(i);
		}
		//System.out.println(Arrays.toString(sequence));
		return sequence;
	}
	
	// setters
	
	// for testing purpose only
	public void setDummySeqInternalNodes (int [] dummySeq) throws Exception{
		for(int i=leafNodeCount; i<nodeCount; i++){
			m_sequences.get(i).setSequence(dummySeq);
		}
	}
	
	// for testing purpose only
	public void setDummyPathBranch(int nodeNr, int seqSite, SubstitutionEvent newEvent){
		m_branches.get(nodeNr).getMutationPath(seqSite).add(newEvent);
	}
	
	// showers ;)
	public void showSequences() {
		Nucleotide nucleo = new Nucleotide();
		for (int i = 0; i < m_nodes.length; i++) {
			int[] tmp_seq = m_sequences.get(i).getSequence();
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
	
	public void showOneSitePath(int seqSite){

		for (int i=0; i<m_branches.size(); i++){

			if(m_branches.get(i).getMutationPath(seqSite).size() == 0){
				System.out.println("no substitution for node:" + i);
			}else{
				for(int j=0; j<m_branches.get(i).getMutationPath(seqSite).size();j++){
					System.out.println("Substituion for node:" + i + " " + m_branches.get(i).getMutationPath(seqSite).get(j).toString());
				}
			}
		}
	}
	
	// getters
	public List<MutableSequence> getSequences(){
		return m_sequences;
	}
	
	public PathBranch getBranch(int branchNr){
		return m_branches.get(branchNr);
	}
	
    /**
     * StateNode implementation
     */
    @Override
    public void setEverythingDirty(boolean bDirty) {
    	super.setEverythingDirty(bDirty);
    	// To-do
    }
    
    /**
     * StateNode implementation 
     * @throws Exception *
     */
    @Override
    protected void store() {
    	super.store();
    	// store branches
    	for(MutableSequence tmpSeq: m_sequences){
    		m_storedsequences.add(tmpSeq.copy());
    	}
    	// store sequences
    	for(PathBranch tmpbranch:m_branches){
    		m_storedbranches.add(tmpbranch.copy());
    	}
    }

    @Override
    public void restore() {
    	super.restore();
    	// restore branches
    	
    	// restore sequences
    }
    
}
