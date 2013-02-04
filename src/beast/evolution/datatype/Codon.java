package beast.evolution.datatype;

import beast.evolution.datatype.DataType.Base;

public class Codon extends Base {

	public Codon () {
		m_nStateCount = 61;
		m_nCodeLength = 3;
		m_sCodeMap = "TTTTTCTTATTGTCTTCCTCATCGTATTACTGTTGCTGGCTTCTCCTACTGCCTCCCCCACCGCATCACCAACAGCGTCGCCGACGGATTATCATAATGACTACCACAACGAATAACAAAAAGAGTAGCAGAAGGGTTGTCGTAGTGGCTGCCGCAGCGGATGACGAAGAGGGTGGCGGAGGG";
		m_mapCodeToStateSet = new int[61][1];
        for (int i = 0; i < 61; i++) {
            m_mapCodeToStateSet[i][0] = i;
        }
	}
	
    @Override
    public String getDescription() {
        return "codon";
    }

}
