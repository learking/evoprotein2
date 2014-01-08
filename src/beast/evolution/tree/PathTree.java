package beast.evolution.tree;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import evoprotein.evolution.datatype.MutableSequence;
import evoprotein.evolution.substitution.SubstitutionEvent;
import beast.core.Description;
import beast.core.Input;
import beast.core.StateNodeInitialiser;
import beast.evolution.alignment.Alignment;
import beast.evolution.datatype.Nucleotide;
import beast.util.TreeParser;


@Description("Tree that enables path sampling")
public class PathTree extends Tree {
	
	public Input<Alignment> m_alignment = new Input<Alignment>("alignment",
			"alignment that contains sequence data");
	
	List<MutableSequence> m_sequences = new ArrayList<MutableSequence>();
	List<MutableSequence> m_storedsequences = new ArrayList<MutableSequence>();
	
	List<PathBranch> m_branches = new ArrayList<PathBranch>();
	List<PathBranch> m_storedbranches = new ArrayList<PathBranch>();
	
	int nucleoSequenceLength;
	int codonSequenceLength;
	
	@Override
	public void initAndValidate() throws Exception {
			// in tree class, this part was, for some reason, commented out
		
			if (m_initial.get() != null && !(this instanceof StateNodeInitialiser)) {
				final Tree other = m_initial.get();
      			root = other.root.copy();
      			nodeCount = other.nodeCount;
      			internalNodeCount = other.internalNodeCount;
      			leafNodeCount = other.leafNodeCount;
      		}
		
			super.initAndValidate();
			
			nucleoSequenceLength = m_alignment.get().getSiteCount();
			codonSequenceLength = nucleoSequenceLength / 3;
			
			for(int i=0; i<nodeCount; i++) {
				m_sequences.add(new MutableSequence(nucleoSequenceLength));
				if(m_nodes[i].isRoot()){
					m_branches.add(new PathBranch());
				}
				else{
					m_branches.add(new PathBranch(nucleoSequenceLength, m_nodes[i].getParent().getNr(), i));
				}
			}
			
			initLeafSequences();
	}
	
    public int getSudoRootNr(){
		int sudoRootNr = 0;
		for (Node childNode : getRoot().getChildren()) {
			if(!childNode.isLeaf()) {
				sudoRootNr = childNode.getNr();
			}
		}
		return sudoRootNr;
    }
	
	private void initLeafSequences() throws Exception{
		for(int i=0; i< leafNodeCount;i++){
			String sequenceID = getNode(i).getID();
			m_sequences.get(i).setSequence(getSequenceByID(sequenceID));
			if(m_sequences.get(i).existStopCodon()){
				throw new Exception("Stop codon encountered in leaf sequence!");
			}
		}
	}
	
	public int[] getSequenceByID(String sequenceID) throws Exception {
		Nucleotide nucleo = new Nucleotide();
		List<Integer> sequenceTarget = new ArrayList<Integer>();
		for (int i = 0; i < m_alignment.get().getNrTaxa(); i++) {
			if (m_alignment.get().sequenceInput.get().get(i).taxonInput.get()
					.toString() == sequenceID) {
				//wrong ! don't try to use member variables
				//sequenceTarget = m_alignment.get().m_pSequences.get().get(i).getSequence(nucleo);
				
				// when we delete insertions, we are operating on m_counts instead of m_pSequences !!
				sequenceTarget = m_alignment.get().getCounts().get(i);
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
			System.out.print("The node is:" + m_nodes[i].getNr() + " ");
			if (m_nodes[i].isLeaf()) {
				System.out.print("The ID is:" + m_nodes[i].getID()
						+ " and length is: " + m_nodes[i].getLength());
			}
			if(!m_nodes[i].isLeaf() && !m_nodes[i].isRoot()){
				System.out.print("Parent:" + m_nodes[i].getParent().getNr() + " Left child:" + m_nodes[i].getLeft().getNr() + " Right child:" + m_nodes[i].getRight().getNr() + " and length is:" + m_nodes[i].getLength());
			}
			if(m_nodes[i].isRoot()){
				System.out.print("This is root. " + "Left child:" + m_nodes[i].getLeft().getNr() + " Right child:" + m_nodes[i].getRight().getNr() + " and length is:" + m_nodes[i].getLength());
			}
			System.out.print(" and the Seq is: " + real_seq + "\n");
		}
	}
	
	public void showPathBranches() throws Exception{
		//show seq path when this branch does not end with sudo-root or root
		for(int i = 0; i < m_branches.size(); i++){
			if(!m_nodes[i].isRoot() && i!=getSudoRootNr()){
				//get seqPath
				int endNodeNr = m_branches.get(i).getEndNodeNr();
				int beginNodeNr = m_branches.get(i).getBeginNodeNr();
				MutableSequence childSeq = m_sequences.get(endNodeNr);
				MutableSequence parentSeq = m_sequences.get(beginNodeNr);
				SeqPath currentSeqPath = m_branches.get(i).getSeqPath(parentSeq, childSeq);
				//print out seqPath so that we can check correctness
				System.out.println("begin node:" + beginNodeNr + " end node:" + endNodeNr);
				System.out.print(currentSeqPath.toString());
			}
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
	
	public List<PathBranch> getBranches(){
		return m_branches;
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
    	// if MCMC works correct, sequences and branches should not be empty
    	// store sequences
    	m_storedsequences.clear();
    	for(MutableSequence tmpSeq: m_sequences){
    		m_storedsequences.add(tmpSeq.copy());
    	}
    	// store branches
    	m_storedbranches.clear();
    	for(PathBranch tmpbranch:m_branches){
    		m_storedbranches.add(tmpbranch.copy());
    	}
    }

    //tmp
    public void fakeStore(){
    	store();
    }
    
    @Override
    public void restore() {
    	super.restore();
    	// restore sequences
    	m_sequences.clear();
    	for(MutableSequence tmpSeq: m_storedsequences){
    		m_sequences.add(tmpSeq.copy());
    	}
    	// restore branches
    	m_branches.clear();
    	for(PathBranch tmpbranch:m_storedbranches){
    		m_branches.add(tmpbranch.copy());
    	}
    	
    	// tmp
    	// showSequences();
    	// System.out.println("#############################################################################");
    }
    
    /**
     * StateNode implementation
     */
    public String toString() {
    	//part I: tree
        String treeStr = super.toString();
        //part II: seqs at nodes
        String seqStr = "";
        for (int i = 0; i < m_sequences.size(); i++) {
        	seqStr += m_sequences.get(i).toString();
        }
        //part III: pathbranchs
        String branchStr = "";
        for (int i = 0; i < m_branches.size(); i++) {
        	//need to implement toString for PathBranch
        	branchStr += m_branches.get(i).toString();
        }
        //Join all parts
        String pathTreeStr = treeStr + "\n" + seqStr + "\n" + branchStr;
        //remove the last character in the string
        pathTreeStr = pathTreeStr.substring(0, pathTreeStr.length() - 1);
        return pathTreeStr;
    }
    
    /**
     * reconstruct PathTree from XML fragment in the form of a DOM node *
     */
    @Override
    public void fromXML(final org.w3c.dom.Node node) {
    	/*************************************************/
    	//get text content
        final String pathTreeStr = node.getTextContent();
        
        /*************************************************/
        //seperate the text into 3 sections: tree, seqs and branches
        String[] parts = pathTreeStr.split("\n");
        String treeStr = parts[0];
        String seqStr = parts[1];
        String branchStr = parts[2];
        
        /*************************************************/
        //parse tree
        //create a tree parser
        //this part needs to be put into a function
        setTreeFromString(treeStr);
        
        /*************************************************/
        //reconstruct m_sequences
        //parse each sequence
        //this part needs to be put into a function too
        setSequencesFromString(seqStr);
        
        /*************************************************/
        //reconstruct m_branches
        setBranchesFromString(branchStr);
        //finish
        
    }

    private void setTreeFromString(String treeStr){
        final TreeParser parser = new TreeParser();
        try {
            parser.thresholdInput.setValue(1e-10, parser);
        } catch (Exception e1) {
            e1.printStackTrace();
        }
        try {
            parser.offsetInput.setValue(0, parser);
            setRoot(parser.parseNewick(treeStr));
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        initArrays();
    }
    
    private void setSequencesFromString(String seqStr){
        List<int[]> allSeqs = new ArrayList<int[]>();
        Pattern seqPattern = Pattern.compile("\\[(.*?)\\]");
        Matcher seqMatcher = seqPattern.matcher(seqStr);
        while (seqMatcher.find()) {
        	//need to remove blank space first
        	String[] tmpStrArr = seqMatcher.group(1).replaceAll(" ", "").split(",");
        	int[] tmpIntArr = new int[tmpStrArr.length];
        	for(int n = 0; n < tmpStrArr.length; n++) {
        	    tmpIntArr[n] = Integer.parseInt(tmpStrArr[n]);
        	}
        	allSeqs.add(tmpIntArr);
        }
        //add seq one by one to m_sequences
        m_sequences.clear();
        for(int i = 0; i < allSeqs.size(); i++){
        	m_sequences.add(new MutableSequence(allSeqs.get(i)));
        }
    }
    
    private void setBranchesFromString(String branchStr){
    	//seperate substrings for diff branches
    	
    	//deal with each branch
    	
    }
    
    /**
     * Loggable implementation
     */    
    public  void init(PrintStream out) throws Exception {
    }
    
    public void log(int nSample, PrintStream out) {
        PathTree thisTree = (PathTree) getCurrent();
        // for each (meaningful) branch
        int rootNr = getRoot().getNr();
        int sudoRootNr = getSudoRootNr();
		for (int branchNr = 0; branchNr < getBranches().size(); branchNr++) {
			if((branchNr != rootNr) && (branchNr != sudoRootNr)) {
				double thisBranchLength = (double) thisTree.getBranch(branchNr).getTotalNumSubstitutions() / (double) nucleoSequenceLength;
				out.print(thisBranchLength + "\t");
			}
		}
		
        int[] dummy = new int[1];
        String sNewick = thisTree.getRoot().toSortedNewick(dummy);
        out.print(sNewick);
        out.print(";");
		
    }
    
    public void close(PrintStream out) {
    }
    
    // tmp
	public void setDummySeqInternalNodesAll (){
		for(int i=leafNodeCount; i<nodeCount; i++){
			int [] newSeq = new int[]{3,3,3,3,3};
			m_sequences.get(i).setSequence(newSeq);
		}
	}
	
	// tmp
	public void setDummySeqInternalNodes (int [] dummySeq) throws Exception{
		for(int i=leafNodeCount; i<nodeCount; i++){
			m_sequences.get(i).setSequence(dummySeq);
		}
	}
	
	// tmp 
	public void setDummySeq(int nodeNr, int [] dummySeq){
		m_sequences.get(nodeNr).setSequence(dummySeq);
	}
	
	// tmp
	public void setDummyPathBranch(int nodeNr, int seqSite, SubstitutionEvent newEvent){
		m_branches.get(nodeNr).getMutationPath(seqSite).add(newEvent);
	}
	
}
