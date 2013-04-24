package beast.evolution.operators;

import beast.core.Description;
import beast.evolution.tree.Node;
import beast.evolution.tree.PathTree;
import beast.util.Randomizer;

@Description("Randomly selects true internal tree node (i.e. not the root or sudo-root) and move node height uniformly in interval " +
        "restricted by the nodes parent and children. At the same time, adjust the substitution times of events within changed branches"+
		"perportionally")
public class BranchLengthOperator extends TreeOperator {

    @Override
    public void initAndValidate() {
    }
	
	@Override
	public double proposal() {
        final PathTree tree = (PathTree) m_tree.get(this);

        // randomly select internal node
        final int nNodeCount = tree.getNodeCount();
        Node node;
        do {
            final int iNodeNr = nNodeCount / 2 + 1 + Randomizer.nextInt(nNodeCount / 2);
            node = tree.getNode(iNodeNr);
            // wrong! also should not be sudo-root
        } while (node.isRoot() || node.isLeaf());
        
        final double fUpper = node.getParent().getHeight();
        final double fLower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
        
        final double newValue = (Randomizer.nextDouble() * (fUpper - fLower)) + fLower;
        // wrong! not only should we set new height for this node, but we should adjust substitution times for events along three branches affected!
        node.setHeight(newValue);
        
        // wrong! should be the complex hastings ratio we discussed!
        return 0.0;
	}

}
