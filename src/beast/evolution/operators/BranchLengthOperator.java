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
        } while (node.isRoot() || node.isLeaf() || node.getHeight()==tree.getRoot().getHeight());
        
        // determine eldest and youngest child 
        int eldestChildNodeNr;
        int youngestChildNodeNr;   
        if(node.getLeft().getHeight() > node.getRight().getHeight()){
        	eldestChildNodeNr = node.getLeft().getNr();
        	youngestChildNodeNr = node.getRight().getNr();
        }else{
        	eldestChildNodeNr = node.getRight().getNr();
        	youngestChildNodeNr = node.getLeft().getNr();        	
        }

        // store node height before change
        final double oldValue = node.getHeight();
        
        // get numbers of substitutions along affected branches
        final int numSubstThisNode = tree.getBranch(node.getNr()).getTotalNumSubstitutions();
        final int numSubstEldest = tree.getBranch(eldestChildNodeNr).getTotalNumSubstitutions();
        final int numSubstYoungest = tree.getBranch(youngestChildNodeNr).getTotalNumSubstitutions();
        
        // get new value
        final double fUpper = node.getParent().getHeight();
        final double fLower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
        final double fYoungest = Math.min(node.getLeft().getHeight(), node.getRight().getHeight());
        final double newValue = (Randomizer.nextDouble() * (fUpper - fLower)) + fLower;
        
        // not only should we set new height for this node, but we should adjust substitution times for events along three branches affected!
        node.setHeight(newValue);
        final double thisNodeScaleFactor = (fUpper - newValue) / (fUpper - oldValue);
        final double eldestChildScaleFactor = (newValue - fLower)/(oldValue - fLower);
        final double youngestChildScaleFactor = (newValue - fYoungest)/(oldValue - fYoungest);
        tree.getBranch(node.getNr()).adjustSubstitutionTimes(thisNodeScaleFactor);
        tree.getBranch(eldestChildNodeNr).adjustSubstitutionTimes(eldestChildScaleFactor);
        tree.getBranch(youngestChildNodeNr).adjustSubstitutionTimes(youngestChildScaleFactor);
        
        /*
        System.out.println("changed node:" + node.getNr() + " num of subst:" + numSubstThisNode + " its eldestChild:" + eldestChildNodeNr + " num of subst:" 
        + numSubstEldest + " its youngest:" + youngestChildNodeNr + " num of subst:" + numSubstYoungest);
        System.out.println("fUpper:" + fUpper + " " + "fLower:" + fLower + " " + "fYoungest:" + fYoungest + " oldValue:" + oldValue + " newValue:" + newValue);
        */
        
        //the complex hastings ratio we discussed
        double logHastingsRatio = (double)numSubstThisNode * (Math.log(fUpper - newValue) - Math.log(fUpper - oldValue)) + 
        		(double)numSubstEldest * (Math.log(newValue - fLower) - Math.log(oldValue - fLower)) + 
        		(double)numSubstYoungest * (Math.log(newValue - fYoungest) - Math.log(oldValue - fYoungest));
        //System.out.println("Hastings ratio:" + logHastingsRatio);
        return logHastingsRatio;
	}

}
