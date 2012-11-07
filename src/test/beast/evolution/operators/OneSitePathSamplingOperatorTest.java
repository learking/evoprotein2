package test.beast.evolution.operators;


import org.junit.Before;
import org.junit.Test;

import test.beast.evoprotein2TestCase;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.operators.OneSitePathSamplingOperator;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.InstantHKY;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.Tree;

import junit.framework.TestCase;

public class OneSitePathSamplingOperatorTest extends evoprotein2TestCase {

	Alignment data;
	Tree tree;
	PathTree pathTree;
	InstantHKY instantHKY;
	RealParameter kappa;
	Frequencies frequences;
	SiteModel siteModel;
	OneSitePathSamplingOperator oneSiteOperator;
	
	@Before
	protected void setUp() throws Exception {
		
		super.setUp();
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
        
        oneSiteOperator = new OneSitePathSamplingOperator();
        oneSiteOperator.initByName("weight", 1.0, "pathtree", pathTree, "siteModel", siteModel);
        
	}

	@Test
	public void testOneSiteOperator() {
		System.out.println((int) 14.99);
	}

}
