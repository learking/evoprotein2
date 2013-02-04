package evoprotein.evolution.datatype;

import java.util.HashSet;

import beast.core.Description;


@Description("used to map seq to codon to codon number")
public class CodonUtil {
	
	public static enum Codon {
		// Stop codons: TAA, TAG, TGA are excluded for now
		
		TTT, TTC, TTA, TTG,
		TCT, TCC, TCA, TCG,
		//TAA, TAG, TAC, TAT,
		TAT, TAC,
		//TGA, TGG, TGC, TGT,
		TGT, TGC,      TGG,
		
		CTT, CTC, CTA, CTG,
		CCT, CCC, CCA, CCG,
		CAT, CAC, CAA, CAG, 
		CGT, CGC, CGA, CGG,
		
		ATT, ATC, ATA, ATG,
		ACT, ACC, ACA, ACG,
		AAT, AAC, AAA, AAG, 
		AGT, AGC, AGA, AGG,
		
		GTT, GTC, GTA, GTG,
		GCT, GCC, GCA, GCG,
		GAT, GAC, GAA, GAG, 
		GGT, GGC, GGA, GGG,
	}
	
	public static enum Nucleotide {
		A , C, G , T
	//  0   1  2   3
	}
	
	/*
	public static int[][][] m_nucleotide2codonMap = new int[][][]{
		{
			{39},{38},{40},{37}
		},
		{
			{23},{},{},{}			
		},
		{
			{},{},{},{}
		},
		{
			{},{},{},{}
		}
	};
	*/
	
	
	
	HashSet<String> codonHashSet;
	
	public CodonUtil(){
		codonHashSet = getCodonHashSet();
	}
	
	public int translate(MutableSequence seq, int startSite) {
		String first = Nucleotide.values()[seq.getSequence()[startSite]].toString();
		String second = Nucleotide.values()[seq.getSequence()[startSite+1]].toString();
		String third = Nucleotide.values()[seq.getSequence()[startSite+2]].toString();
		String thisCodon = first + second + third;
		return Codon.valueOf(thisCodon).ordinal();
	}
	
	int translateOneCodon(int[] codonSeq){
		return -1;
	}
	
	public static HashSet<String> getCodonHashSet() {

		  HashSet<String> codonHashSet = new HashSet<String>();

		  for (Codon c : Codon.values()) {
		      codonHashSet.add(c.name());
		  }

		  return codonHashSet;
	}
	
	public String int2Codon(int firstNucleotide, int secondNucleotide, int thirdNucleotide){
		String first = Nucleotide.values()[firstNucleotide].toString();
		String second = Nucleotide.values()[secondNucleotide].toString();
		String third = Nucleotide.values()[thirdNucleotide].toString();
		String thisCodon = first + second + third;
		return thisCodon;
	}
	
	public boolean containsCodon(String codon){
		return codonHashSet.contains(codon);
	}
	
	
}
