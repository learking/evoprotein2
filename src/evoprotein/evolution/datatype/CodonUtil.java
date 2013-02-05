package evoprotein.evolution.datatype;

import java.util.Arrays;
import java.util.HashMap;
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
	

	int[][][] m_nucleotide2codonMap;
	
	//public static HashMap<int[],Integer> m_nucleotide2codonMap;
	
	HashSet<String> codonHashSet;
	
	public CodonUtil(){
		codonHashSet = getCodonHashSet();
		m_nucleotide2codonMap = createNucleotide2CodonMap();
	}
	
	int[][][] createNucleotide2CodonMap(){
		int[][][] nucleotide2codonMap = new int[4][4][4];
		nucleotide2codonMap[3][3][3]=0;
		nucleotide2codonMap[3][3][1]=1;
		nucleotide2codonMap[3][3][0]=2;
		nucleotide2codonMap[3][3][2]=3;
		
		nucleotide2codonMap[3][1][3]=4;
		nucleotide2codonMap[3][1][1]=5;
		nucleotide2codonMap[3][1][0]=6;
		nucleotide2codonMap[3][1][2]=7;
		
		nucleotide2codonMap[3][0][3]=8;
		nucleotide2codonMap[3][0][1]=9;
		nucleotide2codonMap[3][0][0]=-1; //stop codon
		nucleotide2codonMap[3][0][2]=-1; //stop codon
		
		nucleotide2codonMap[3][2][3]=10;
		nucleotide2codonMap[3][2][1]=11;
		nucleotide2codonMap[3][2][0]=-1; //stop codon
		nucleotide2codonMap[3][2][2]=12;
		
		nucleotide2codonMap[1][3][3]=13;
		nucleotide2codonMap[1][3][1]=14;
		nucleotide2codonMap[1][3][0]=15;
		nucleotide2codonMap[1][3][2]=16;
		
		nucleotide2codonMap[1][1][3]=17;
		nucleotide2codonMap[1][1][1]=18;
		nucleotide2codonMap[1][1][0]=19;
		nucleotide2codonMap[1][1][2]=20;
		
		nucleotide2codonMap[1][0][3]=21;
		nucleotide2codonMap[1][0][1]=22;
		nucleotide2codonMap[1][0][0]=23;
		nucleotide2codonMap[1][0][2]=24;
		
		nucleotide2codonMap[1][2][3]=25;
		nucleotide2codonMap[1][2][1]=26;
		nucleotide2codonMap[1][2][0]=27;
		nucleotide2codonMap[1][2][2]=28;
		
		nucleotide2codonMap[0][3][3]=29;
		nucleotide2codonMap[0][3][1]=30;
		nucleotide2codonMap[0][3][0]=31;
		nucleotide2codonMap[0][3][2]=32;

		nucleotide2codonMap[0][1][3]=33;
		nucleotide2codonMap[0][1][1]=34;
		nucleotide2codonMap[0][1][0]=35;
		nucleotide2codonMap[0][1][2]=36;
		
		nucleotide2codonMap[0][0][3]=37;
		nucleotide2codonMap[0][0][1]=38;
		nucleotide2codonMap[0][0][0]=39;
		nucleotide2codonMap[0][0][2]=40;
		
		nucleotide2codonMap[0][2][3]=41;
		nucleotide2codonMap[0][2][1]=42;
		nucleotide2codonMap[0][2][0]=43;
		nucleotide2codonMap[0][2][2]=44;
		
		nucleotide2codonMap[2][3][3]=45;
		nucleotide2codonMap[2][3][1]=46;
		nucleotide2codonMap[2][3][0]=47;
		nucleotide2codonMap[2][3][2]=48;
		
		nucleotide2codonMap[2][1][3]=49;
		nucleotide2codonMap[2][1][1]=50;
		nucleotide2codonMap[2][1][0]=51;
		nucleotide2codonMap[2][1][2]=52;
		
		nucleotide2codonMap[2][0][3]=53;
		nucleotide2codonMap[2][0][1]=54;
		nucleotide2codonMap[2][0][0]=55;
		nucleotide2codonMap[2][0][2]=56;
		
		nucleotide2codonMap[2][2][3]=57;
		nucleotide2codonMap[2][2][1]=58;
		nucleotide2codonMap[2][2][0]=59;
		nucleotide2codonMap[2][2][2]=60;
		return nucleotide2codonMap;
	}
	
	
	public int translate(MutableSequence seq, int startSite){
		int[] codonSeq = Arrays.copyOfRange(seq.getSequence(), startSite, startSite + 3);
		return m_nucleotide2codonMap[codonSeq[0]][codonSeq[1]][codonSeq[2]];
	}
	
	
	/*
	public int translate(MutableSequence seq, int startSite) {
		String first = Nucleotide.values()[seq.getSequence()[startSite]].toString();
		String second = Nucleotide.values()[seq.getSequence()[startSite+1]].toString();
		String third = Nucleotide.values()[seq.getSequence()[startSite+2]].toString();
		String thisCodon = first + second + third;
		return Codon.valueOf(thisCodon).ordinal();
	}
	*/
	
	public static HashSet<String> getCodonHashSet() {

		  HashSet<String> codonHashSet = new HashSet<String>();

		  for (Codon c : Codon.values()) {
		      codonHashSet.add(c.name());
		  }

		  return codonHashSet;
	}
	
	/*
	public String int2Codon(int firstNucleotide, int secondNucleotide, int thirdNucleotide){
		String first = Nucleotide.values()[firstNucleotide].toString();
		String second = Nucleotide.values()[secondNucleotide].toString();
		String third = Nucleotide.values()[thirdNucleotide].toString();
		String thisCodon = first + second + third;
		return thisCodon;
	}
	*/
	
	public boolean containsCodon(String codon){
		return codonHashSet.contains(codon);
	}
	
	
}
