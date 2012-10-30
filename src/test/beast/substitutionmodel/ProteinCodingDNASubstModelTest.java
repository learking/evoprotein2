package test.beast.substitutionmodel;


import org.junit.Before;
import org.junit.Test;

import evoprotein.evolution.datatype.MutableSequence;

import beast.evolution.substitutionmodel.ProteinCodingDNASubstModel;

import junit.framework.TestCase;

public class ProteinCodingDNASubstModelTest extends TestCase {

	ProteinCodingDNASubstModel ourModel;
	MutableSequence seq;
	
	@Before
	protected void setUp() throws Exception {
		super.setUp();
		
		ourModel = new ProteinCodingDNASubstModel();
		
		seq = new MutableSequence(3);
		int [] intSeq = new int [] {1,3,0};
		seq.setSequence(intSeq);
	}

	@Test
	public void testGetSubstAwayRate() throws Exception {
		ourModel.getSubstAwayRate(seq);
	}

}
