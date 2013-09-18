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
        if (userDataTypeInput.get() != null) {
            m_dataType = userDataTypeInput.get();
        } else {
            if (types.indexOf(dataTypeInput.get()) < 0) {
                throw new Exception("data type + '" + dataTypeInput.get() + "' cannot be found. " +
                        "Choose one of " + Arrays.toString(types.toArray(new String[0])));
            }
            List<String> sDataTypes = AddOnManager.find(beast.evolution.datatype.DataType.class, IMPLEMENTATION_DIR);
            for (String sDataType : sDataTypes) {
                DataType dataType = (DataType) Class.forName(sDataType).newInstance();
                if (dataTypeInput.get().equals(dataType.getDescription())) {
                    m_dataType = dataType;
                    break;
                }
            }
        }

        // grab data from child sequences
        taxaNames.clear();
        stateCounts.clear();
        counts.clear();
        for (Sequence seq : sequenceInput.get()) {
            //m_counts.add(seq.getSequence(getMap()));
            counts.add(seq.getSequence(m_dataType));
            if (taxaNames.indexOf(seq.taxonInput.get()) >= 0) {
                throw new Exception("Duplicate taxon found in alignment: " + seq.taxonInput.get());
            }
            taxaNames.add(seq.taxonInput.get());
            stateCounts.add(seq.totalCountInput.get());
        }
        if (counts.size() == 0) {
            // no sequence data
            throw new Exception("Sequence data expected, but none found");
        }

        // handle insertions first, because position numbering will change after we deletion things
        removeInsertions();
        
        // mark down positions where deletions happen
        if(noInsertion(counts.get(0))){
            // store all deletion-caused gap positions first
        	findDeletions();
        	// remove all deletion-caused gap positions
        	removeDeletions();
        }else{
        	throw new Exception("there shouldn't be any insertions left!");
        }
        
        // Sanity check: make sure sequences are of same length
        int nLength = counts.get(0).size();
        for (List<Integer> seq : counts) {
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
    	for(int nTaxa = 0; nTaxa < counts.size(); nTaxa++){
    		newSequences.add(new ArrayList<Integer>());
    	}
    	
    	List<Integer> referenceSeq = counts.get(0);
    	// if gap in the first (reference) seq, delete same positions in other seqs
    	for (int AAsite=0; AAsite < referenceSeq.size(); AAsite += 3) {
    		if(referenceSeq.get(AAsite) != GAP_INT){
    			addAA(newSequences, counts, AAsite);
    		}
    	}
    	counts = newSequences;
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
    	for(int i = 0 ; i < counts.get(0).size(); i +=3){
    		if(isDeletion(counts, i)){
    			// add i to the collection of deletion positions
    			m_deletionPositions.add(i);
    		}
    	}
    }
    
    void removeDeletions() {
    	
    	if(!m_deletionPositions.isEmpty()){
    		
        	List<List<Integer>> newSequences = new ArrayList<List<Integer>>();
        	for(int nTaxa = 0; nTaxa < counts.size(); nTaxa++){
        		newSequences.add(new ArrayList<Integer>());
        	}
        	
        	// if this position is not in m_deletionPositions, add it to newSequences
        	List<Integer> referenceSeq = counts.get(0);
        	// if gap in the first (reference) seq, delete same positions in other seqs
        	for (int AAsite=0; AAsite < referenceSeq.size(); AAsite += 3) {
        		// if not a deletion site
        		if(!m_deletionPositions.contains(AAsite)){
        			addAA(newSequences, counts, AAsite);
        		}
        	}
        	counts = newSequences;
        	
    	}
    }
    
    // getter for m_deletionPositions
    public Set<Integer> getDeletionPositions(){
    	return m_deletionPositions;
    }
    
}


