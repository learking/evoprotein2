package test.beast.evolution.operators;


import java.util.List;

import org.junit.Test;

import evoprotein.evolution.substitution.SubstitutionEvent;

import test.beast.evoprotein2TestCase;
import test.beast.core.testMCMC;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.PathTreeLikelihood;
import beast.evolution.operators.PathSamplingOperator;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.InstantHKY;
import beast.evolution.tree.Node;
import beast.evolution.tree.PathBranch;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.Tree;

public class PathSamplingOperatorTest extends evoprotein2TestCase {
	
	Alignment data;
	Tree tree;
	PathTree pathTree;
	InstantHKY instantHKY;
	RealParameter kappa;
	Frequencies frequences;
	SiteModel siteModel;
	PathSamplingOperator pathSamplingOperator;
	
	@Test
	public void testPathSamplingOperator() throws Exception {
        data = getDummyAlignment();
        //data = getAlignment();
		
        tree = getTree(data);
		
        pathTree = new PathTree();
		pathTree.initByName("initial", tree, "alignment", data);
        
		kappa = getKappa("1.5");
        
        frequences = new Frequencies();
        frequences.initByName("data", data);
        
        instantHKY = new InstantHKY();
        instantHKY.initByName("kappa", kappa, "frequencies", frequences);
        
        siteModel = new SiteModel();
        siteModel.initByName("gammaCategoryCount", 1, "substModel", instantHKY);
        
        pathSamplingOperator = new PathSamplingOperator();
        pathSamplingOperator.initByName("weight", 1.0, "pathtree", pathTree, "siteModel", siteModel);
        
        
		int rootNr = pathTree.getRoot().getNr();
		int sudoRootNr = 0;
		for (Node childNode : pathTree.getRoot().getChildren()) {
			if(!childNode.isLeaf()) {
				sudoRootNr = childNode.getNr();
			}
		}
        
		// simulate multiple runs inside MCMC
		for (int iSample = 0; iSample < 3 ; iSample++){
			// Pupko
			for (int seqSite = 0; seqSite < pathTree.getSequences().get(0).getSequence().length; seqSite ++) {
				pathSamplingOperator.PupkoOneSite(pathTree, seqSite);
			}
			// pathSamplingOperator.setPathLogDenstiyToZero();
			// Nielsen
			for (int branchNr = 0; branchNr < pathTree.getBranches().size(); branchNr++) {
				if((branchNr != rootNr) && (branchNr != sudoRootNr)) {
					pathSamplingOperator.NielsenSampleOneBranch(pathTree, branchNr);
				}
			}
			System.err.println("HastingRatio:" + pathSamplingOperator.getPathLogDensity());
			
			// function here to check the correctness of the above two methods:
			checkPathSampling(pathTree, sudoRootNr, rootNr);
			pathSamplingOperator.setPathLogDenstiyToZero();
			System.out.println("===============================================");
		}
		
	}

	public void checkPathSampling(PathTree pathTree, int sudoRootNr, int rootNr) throws Exception{
		// for each branch (except sudoRoot and root ended branches)
		for (int branchNr = 0; branchNr < pathTree.getBranches().size() ; branchNr++) {
			if ((branchNr != sudoRootNr) && (branchNr != rootNr)){
				// for each site within branch
				PathBranch currentBranch = pathTree.getBranch(branchNr);
				int beginNodeNr = currentBranch.getBeginNodeNr();
				int endNodeNr = currentBranch.getEndNodeNr();
				int [] beginSeq = pathTree.getSequences().get(beginNodeNr).getSequence();
				int [] endSeq = pathTree.getSequences().get(endNodeNr).getSequence();
				for (int seqSite = 0; seqSite < pathTree.getSequences().get(0).getSequence().length; seqSite++){
					// check whether begin-end fit or last-nucleo-end fit
					List<SubstitutionEvent> eventsAtThisSite = currentBranch.getMutationPath(seqSite);
					int beginNucleotide = beginSeq[seqSite];
					int endNucleotide = endSeq[seqSite];
					if(eventsAtThisSite.size() == 0){
						if(beginNucleotide != endNucleotide){
							throw new Exception("begin != end");
						}
					}
					else{
						int lastNucleotide = eventsAtThisSite.get(eventsAtThisSite.size() - 1).getCurrentNucleotide();
						if(lastNucleotide != endNucleotide){
							System.out.println("last:" + lastNucleotide + " end:" + endNucleotide + 
									" begin:" + beginNucleotide + " branchNr:" + branchNr + " seqSite:" + seqSite
									+ "events:");
							for (SubstitutionEvent tmpEvent : eventsAtThisSite){
								System.out.print(tmpEvent.toString());
							}
							throw new Exception("last != end");
						}
					}
				}
			}
		}
		
	}
	
}
