package evoprotein.evolution.datatype;

import java.util.HashSet;

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
		
		// Stop codons: TAA, TAG, TGA are excluded for now
		
		//TAA, TAG, TAC, TAT,
		TAC, TAT,
		//TGA, TGG, TGC, TGT,
		TGG, TGC, TGT,
		TCA, TCG, TCC, TCT,
		TTA, TTG, TTC, TTT,
	}
	
	public enum Nucleotide {
		A , C, G , T
	//  0   1  2   3
	}
	
	HashSet<String> codonHashSet;
	
	public CodonUtil(){
		codonHashSet = getCodonHashSet();
	}
	
	public int translate(MutableSequence seq, int startSite) throws Exception {
		String first = Nucleotide.values()[seq.getSequence()[startSite]].toString();
		String second = Nucleotide.values()[seq.getSequence()[startSite+1]].toString();
		String third = Nucleotide.values()[seq.getSequence()[startSite+2]].toString();
		String thisCodon = first + second + third;
		if(containsCodon(thisCodon)){
			return Codon.valueOf(thisCodon).ordinal();
		}else{
			throw new Exception("Encounter a stop codon!");
		}
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
