package test.beast.evolution.operators;


import org.junit.Test;

import test.beast.evoprotein2TestCase;

import beast.evolution.alignment.Alignment;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.Tree;

import junit.framework.TestCase;

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
	public void test() {
		fail("Not yet implemented");
	}

}
