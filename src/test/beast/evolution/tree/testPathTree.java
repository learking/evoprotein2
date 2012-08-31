package test.beast.evolution.tree;


import org.junit.Test;

import test.beast.evoprotein2TestCase;

import beast.core.Description;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.Tree;

@Description("Test PathTree")
public class testPathTree extends evoprotein2TestCase {

	Alignment data;
	Tree tree;
	
    @Override
    protected void setUp() throws Exception {
        super.setUp();
        data = getAlignment();
        tree = getTree(data);
    }
	
	@Test
	public void testTreeInitialization() throws Exception{
		PathTree pathTree = new PathTree();
		pathTree.initByName("initial", tree, "alignment", data);
		int[][] sequences = pathTree.getSequences();
		assertEquals(sequences[0].length, 3000);
	}

}
