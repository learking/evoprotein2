/**
 * 
 */
package beast.evolution.alignment;

import java.util.Arrays;
import java.util.List;

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

        // handle insertions and deletions here
        findInDels();
        removeInsertions();
        
        // Sanity check: make sure sequences are of same length
        int nLength = m_counts.get(0).size();
        for (List<Integer> seq : m_counts) {
            if (seq.size() != nLength) {
                throw new Exception("Two sequences with different length found: " + nLength + " != " + seq.size());
            }
        }
        
        calcPatterns();
    }
	
    void findInDels(){
    	// how to find indel at the beginning, middle section and end region?
    	
    }
    
    void removeInsertions(){
    	// do something to sequences stored in m_counts to remove insertions
    	
    }

}
