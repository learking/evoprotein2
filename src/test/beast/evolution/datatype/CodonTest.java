package test.beast.evolution.datatype;

import static org.junit.Assert.*;

import org.junit.Test;

import beast.evolution.datatype.Codon;

public class CodonTest {

	@Test
	public void testString2State() {
		Codon codon = new Codon();
		int [] stateSeq = new int[]{60,0,1};
		System.out.println(codon.state2string(stateSeq));
	}

}
