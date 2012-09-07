package test.beast.evolution.operators;


import org.junit.Test;

import test.beast.evoprotein2TestCase;

import beast.evolution.alignment.Alignment;
import beast.evolution.operators.PathSamplingOperator;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.Tree;

public class PathSamplingOperatorTest extends evoprotein2TestCase {

	Alignment data;
	Tree tree;
	PathTree pathTree;
	
    @Override
    protected void setUp() throws Exception {
        super.setUp();
        data = getAlignment();
        tree = getTree(data);
		PathTree pathTree = new PathTree();
		pathTree.initByName("initial", tree, "alignment", data);
    }
	
	@Test
	public void testPathSamplingOperator() throws Exception {
		PathSamplingOperator pathSamplingOperator = new PathSamplingOperator();
		pathSamplingOperator.initByName("pathtree", pathTree);
	}

}
