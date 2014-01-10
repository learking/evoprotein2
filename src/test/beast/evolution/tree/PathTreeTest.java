package test.beast.evolution.tree;


import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.util.List;

import javax.xml.parsers.DocumentBuilderFactory;

import junit.framework.TestCase;

import org.junit.Ignore;
import org.junit.Test;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import evoprotein.evolution.datatype.MutableSequence;
import evoprotein.evolution.substitution.SubstitutionEvent;
import test.beast.evoprotein2TestCase;
import beast.core.Description;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;

@Description("Test PathTree")
public class PathTreeTest extends TestCase {

	Alignment data;
	Tree tree;
	
    static public Alignment getSimpleAlignment() throws Exception {
        Sequence Human = new Sequence("Human", "CCA");
        Sequence Chimpanzee = new Sequence("Chimpanzee", "CCC");
        Sequence Gibbon = new Sequence("Gibbon", "CTA");
        Sequence Gorilla = new Sequence("Gorilla", "CCA");

        Alignment data = new Alignment();
        data.initByName("sequence", Human, "sequence", Chimpanzee, "sequence", Gorilla, "sequence", Gibbon,
                "dataType", "nucleotide"
        );
        return data;
    }
	
    static public Tree getSimpleTree(Alignment data) throws Exception {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", "(((Human: 0.1, Chimpanzee: 0.1): 0.4, Gibbon:0.5):0 , Gorilla: 0.5)");
        return tree;
    }
	
    @Override
    protected void setUp() throws Exception {
    	super.setUp();
        data = getSimpleAlignment();
        tree = getSimpleTree(data);
    }
	
    @Test
	@Ignore public void testToString() throws Exception{
		PathTree pathTree = new PathTree();
		pathTree.initByName("initial", tree, "alignment", data);
		//initialize internal node seq
		pathTree.setDummySeq(4, new int[]{1,1,1});
		pathTree.setDummySeq(5, new int[]{1,0,1});
		pathTree.setDummySeq(6, new int[]{1,0,1});
		//add mutationpaths for each branch
		pathTree.showSequences();
		//add branches
		//path ends at node 0, at site 2 (0,1,2), we add a substitution "begin C(1), end A(0), interval 0.03" 
		pathTree.setDummyPathBranch(0, 2, new SubstitutionEvent(1,0,0.03));
		
		pathTree.setDummyPathBranch(1, 2, new SubstitutionEvent(1,0,0.02));
		pathTree.setDummyPathBranch(1, 2, new SubstitutionEvent(0,1,0.04));
		
		pathTree.setDummyPathBranch(3, 1, new SubstitutionEvent(0,3,0.1));
		pathTree.setDummyPathBranch(3, 2, new SubstitutionEvent(1,0,0.3));
		
		pathTree.setDummyPathBranch(2, 2, new SubstitutionEvent(1,0,0.2));
		pathTree.setDummyPathBranch(2, 1, new SubstitutionEvent(0,1,0.4));
		
		//don't forget to add subs between node 4 and 5
		pathTree.setDummyPathBranch(4, 1, new SubstitutionEvent(0,1,0.05));
		
		//show all pathBranches to verify that the pathTree is what we designed
		pathTree.showPathBranches();
		System.out.print(pathTree.toString());
	}
	
    @Test
    public void testFromXML() throws Exception{
    	//prepare parser
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        
        //under linux
        //String stateFileName = "/home/kuangyu/Dropbox/research/meeting_memo/Jan_2014/pathTree.xml.state";
        //under mac
        String stateFileName = "/Users/kwang2/Desktop/Dropbox/research/meeting_memo/Jan_2014/pathTree.xml.state";
        
        //String stateFileName = "/home/kuangyu/workspace/beast2/realData_twoStruct_gapped.xml.state";
        Document doc = factory.newDocumentBuilder().parse(new File(stateFileName));
        doc.normalize();
        //get node that belongs to PathTree
        final NodeList nodes = doc.getElementsByTagName("*");
        final Node topNode = nodes.item(0);
        final NodeList children = topNode.getChildNodes();  
        final Node child = children.item(1);

        //verify its ID:
        //final String sID = child.getAttributes().getNamedItem("id").getNodeValue();
        //System.out.println(sID);

        //show its text content:
        PathTree pathTree = new PathTree();
        pathTree.fromXML(child);

        //
		pathTree.showPathBranches();
		System.out.print(pathTree.toString());
    }
    
}
