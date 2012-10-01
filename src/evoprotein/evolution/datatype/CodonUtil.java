package evoprotein.evolution.datatype;

import beast.core.Description;


@Description("used to map seq to codon to codon number")
public class CodonUtil {
	
	public static enum Codon {
		AAA, AAG, AAC, AAT, 
		AGA, AGG, AGC, AGT,
		ACA, ACG, ACC, ACT,
		ATA, ATG, ATC, ATT,
		
		GAA, GAG, GAC, GAT, 
		GGA, GGG, GGC, GGT,
		GCA, GCG, GCC, GCT,
		GTA, GTG, GTC, GTT,
		
		CAA, CAG, CAC, CAT, 
		CGA, CGG, CGC, CGT,
		CCA, CCG, CCC, CCT,
		CTA, CTG, CTC, CTT,
		
		TAA, TAG, TAC, TAT, 
		TGA, TGG, TGC, TGT,
		TCA, TCG, TCC, TCT,
		TTA, TTG, TTC, TTT,
	}
	
	public enum Nucleotide {
		A , C, G , T
	//  0   1  2   3
	}
	
	public int translate(MutableSequence seq, int startSite) {
		String first = Nucleotide.values()[seq.getSequence()[startSite]].toString();
		String second = Nucleotide.values()[seq.getSequence()[startSite+1]].toString();
		String third = Nucleotide.values()[seq.getSequence()[startSite+2]].toString();
		String thisCodon = first + second + third;
		return Codon.valueOf(thisCodon).ordinal();
	}
	
}
