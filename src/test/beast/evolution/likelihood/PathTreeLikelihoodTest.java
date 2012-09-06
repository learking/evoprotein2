/**
 * 
 */
package test.beast.evolution.likelihood;


import org.junit.Test;

import evoprotein.evolution.substitution.SubstitutionEvent;

import test.beast.BEASTTestCase;
import test.beast.evoprotein2TestCase;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.PathTreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.InstantHKY;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.Tree;

/**
 * @author kuangyu
 *
 */
public class PathTreeLikelihoodTest extends evoprotein2TestCase {

	// weekend task 1:
	// create a path tree with internal states assigned and path mapped out
	// use this tree to test pathTreeLikelihood
	
	Alignment data;
	Tree tree;
	
    @Override
    protected void setUp() throws Exception {
        super.setUp();
        data = getDummyAlignment();
        tree = getTree(data);
    }
	
	@Test
	public void testLikelihoodCalculation() throws Exception {
		// data has been setup already in setup()
		
		// set up PathTree
		
		PathTree pathTree = new PathTree();
		pathTree.initByName("initial", tree, "alignment", data);
		
		/*
		int nodeNr = 5;
		int leftNr = pathTree.getNode(nodeNr).getLeft().getNr();
		int rightNr = pathTree.getNode(nodeNr).getRight().getNr();
		System.out.println("left nr:" + leftNr + " LH:"+ pathTree.getNode(nodeNr).getLeft().getHeight() + " right Nr:" + rightNr + " RH:"+ pathTree.getNode(nodeNr).getRight().getHeight() + " height is:" + pathTree.getNode(nodeNr).getHeight());
		*/
		
		int [] dummySeq = {1,1,0,0,0};
		pathTree.setDummySeqInternalNodes(dummySeq);
		pathTree.showSequences();
		pathTree.setDummyPathBranch(2, 1, new SubstitutionEvent(1,3,0.25));
		pathTree.setDummyPathBranch(3, 1, new SubstitutionEvent(1,3,0.2));
		pathTree.setDummyPathBranch(4, 1, new SubstitutionEvent(1,3,0.2));
		
		//pathTree.showOneSitePath(1);
		
		// set up Substitution model
		RealParameter f = new RealParameter(new Double[]{0.25, 0.25, 0.25, 0.25});
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", f, "estimate", false);
    
        InstantHKY instantHKY = new InstantHKY();
        instantHKY.initByName("kappa", "2.0", "frequencies", freqs);
		
		// set up SiteModel (which includes setting up substitution model)
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", instantHKY);
		
		// create PathTreeLikelihood
        PathTreeLikelihood likelihood = new PathTreeLikelihood();
        likelihood.initByName("PathTree", pathTree, "siteModel", siteModel);
		
		// test PathTreeLikelihood's correctness
        double oneSiteLogP = 0;
        oneSiteLogP = likelihood.calculateOneSiteLogP(1);
        
        // Tasks for today, Sep 5th, 2012
        // 2. implement oneSiteLikelihood for "PathTreeLikelihood" and test it
        // 3. extend oneSiteLikelihood to all sites
        
        assertEquals(oneSiteLogP, 4.0, BEASTTestCase.PRECISION);
        //assertEquals(oneSiteLogP, -3.189441542, BEASTTestCase.PRECISION);
        
	}

}
