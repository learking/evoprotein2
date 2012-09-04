package test.beast.evolution.tree;


import java.util.List;

import org.junit.Test;

import evoprotein.evolution.datatype.MutableSequence;

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
		List<MutableSequence> sequences = pathTree.getSequences();
		assertEquals(sequences.get(0).getSequence().length, 3000);
	}

}
