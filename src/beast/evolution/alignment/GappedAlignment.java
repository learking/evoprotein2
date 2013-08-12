/**
 * 
 */
package beast.evolution.alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beast.evolution.datatype.DataType;
import beast.util.AddOnManager;

/**
 * @author kuangyu
 *
 */
public class GappedAlignment extends Alignment {

	/**
	 * 
	 */
	public GappedAlignment() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * constructor for testing purpose still
	 * @param sequences
	 * @param stateCount
	 * @param dataType
	 * @throws Exception
	 */
	public GappedAlignment(List<Sequence> sequences, Integer stateCount,
			String dataType) throws Exception {
		super(sequences, stateCount, dataType);
		// TODO Auto-generated constructor stub
	}

    @Override
    public void initAndValidate() throws Exception {
        // determine data type, either user defined or one of the standard ones
        if (m_userDataType.get() != null) {
            m_dataType = m_userDataType.get();
        } else {
            if (m_sTypes.indexOf(m_sDataType.get()) < 0) {
                throw new Exception("data type + '" + m_sDataType.get() + "' cannot be found. " +
                        "Choose one of " + Arrays.toString(m_sTypes.toArray(new String[0])));
            }
            List<String> sDataTypes = AddOnManager.find(beast.evolution.datatype.DataType.class, IMPLEMENTATION_DIR);
            for (String sDataType : sDataTypes) {
                DataType dataType = (DataType) Class.forName(sDataType).newInstance();
                if (m_sDataType.get().equals(dataType.getDescription())) {
                    m_dataType = dataType;
                    break;
                }
            }
        }

        // grab data from child sequences
        m_sTaxaNames.clear();
        m_nStateCounts.clear();
        m_counts.clear();
        for (Sequence seq : m_pSequences.get()) {
            //m_counts.add(seq.getSequence(getMap()));
            m_counts.add(seq.getSequence(m_dataType));
            if (m_sTaxaNames.indexOf(seq.m_sTaxon.get()) >= 0) {
                throw new Exception("Duplicate taxon found in alignment: " + seq.m_sTaxon.get());
            }
            m_sTaxaNames.add(seq.m_sTaxon.get());
            m_nStateCounts.add(seq.m_nTotalCount.get());
        }
        if (m_counts.size() == 0) {
            // no sequence data
            throw new Exception("Sequence data expected, but none found");
        }

        // handle insertions first, because position numbering will change after we deletion things
        removeInsertions();
        
        // mark down positions where deletions happen
        if(noInsertion(m_counts.get(0))){
            findDeletions();
        }else{
        	throw new Exception("there shouldn't be any insertions left!");
        }
        
        // Sanity check: make sure sequences are of same length
        int nLength = m_counts.get(0).size();
        for (List<Integer> seq : m_counts) {
            if (seq.size() != nLength) {
                throw new Exception("Two sequences with different length found: " + nLength + " != " + seq.size());
            }
        }
        
        calcPatterns();
    }
    
    // handle indels
    final static public int GAP_INT = 16;
    
    boolean noInsertion(List<Integer> referenceSeq){
    	boolean noInsertionFlag = true;
    	for(int i = 0; i < referenceSeq.size(); i++){
    		if(referenceSeq.get(i) == GAP_INT){
    			noInsertionFlag = false;
    			break;
    		}
    	}
    	return noInsertionFlag;
    }

    void removeInsertions(){
    	List<List<Integer>> newSequences = new ArrayList<List<Integer>>();
    	for(int nTaxa = 0; nTaxa < m_counts.size(); nTaxa++){
    		newSequences.add(new ArrayList<Integer>());
    	}
    	
    	List<Integer> referenceSeq = m_counts.get(0);
    	// if gap in the first (reference) seq, delete same positions in other seqs
    	for (int AAsite=0; AAsite < referenceSeq.size(); AAsite += 3) {
    		if(referenceSeq.get(AAsite) != GAP_INT){
    			addAA(newSequences, m_counts, AAsite);
    		}
    	}
    	m_counts = newSequences;
    }
    
    void addAA(List<List<Integer>> newSequences, List<List<Integer>> oldSequences, int AAsite){
    	for (int nTaxa = 0; nTaxa < oldSequences.size(); nTaxa++) {
    		newSequences.get(nTaxa).add(oldSequences.get(nTaxa).get(AAsite));
    		newSequences.get(nTaxa).add(oldSequences.get(nTaxa).get(AAsite + 1));
    		newSequences.get(nTaxa).add(oldSequences.get(nTaxa).get(AAsite + 2));    		
    	}
    }
    
    boolean isDeletion(List<List<Integer>> sequences, int AAsite) throws Exception{
    	boolean deletionFlag = false;
    	for(int nTaxa = 1; nTaxa < sequences.size(); nTaxa++){
    		List<Integer> tmpSeq = sequences.get(nTaxa);
    		if(tmpSeq.get(AAsite) == GAP_INT && tmpSeq.get(AAsite + 1) == GAP_INT && tmpSeq.get(AAsite + 2) == GAP_INT){
    			deletionFlag = true;
    			break;
    		}else if (tmpSeq.get(AAsite) != GAP_INT && tmpSeq.get(AAsite + 1) != GAP_INT && tmpSeq.get(AAsite + 2) != GAP_INT){
    			continue;
    		}else{
    			throw new Exception("A deletion should be 3 nucleotide wide only");
    		}
    	}
    	return deletionFlag;
    }
    
    //shouldn't be a list, but a set
    protected Set<Integer> m_deletionPositions = new HashSet<Integer>();
    
	// if gap in other seq/seqs mark down the positions
    void findDeletions() throws Exception{
    	// skip the reference seq
    	for(int i = 0 ; i < m_counts.get(0).size(); i +=3){
    		if(isDeletion(m_counts, i)){
    			// add i to the collection of deletion positions
    			m_deletionPositions.add(i);
    		}
    	}
    }
    
    // getter for m_deletionPositions
    public Set<Integer> getDeletionPositions(){
    	return m_deletionPositions;
    }
    
}


