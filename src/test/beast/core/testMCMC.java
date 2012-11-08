/**
 * 
 */
package test.beast.core;


import java.util.ArrayList;

import org.junit.Test;

import evoprotein.proteinstructure.SolventAccessibility;

import beast.core.Distribution;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Plugin;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.PathTreeLikelihood;
import beast.evolution.operators.AllSitesPathSamplingOperator;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.InstantHKY;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.Tree;

import test.beast.evoprotein2TestCase;

/**
 * @author kuangyu
 *
 */
public class testMCMC extends evoprotein2TestCase {

	Alignment data;
	Tree tree;
	PathTree pathTree;
	InstantHKY instantHKY;
	InstantHKY proposalInstantHKY;
	RealParameter kappa;
	RealParameter proposalKappa;
	Frequencies frequencies;
	Frequencies proposalFrequencies;
	SiteModel siteModel;
	SiteModel proposalSiteModel;
	PathTreeLikelihood pathTreeLikelihood;
	Distribution likelihood;

	Logger tmpLogger;
	ArrayList<Plugin> log;
	
	AllSitesPathSamplingOperator allSitesPathSamplingOperator;
	
	//SolventAccessibility solventAccessibility;
	
    @Override
    protected void setUp() throws Exception {
        super.setUp();
        
        data = getAlignment();
        
        tree = getTree(data);
		
        pathTree = new PathTree();
		pathTree.initByName("initial", tree, "alignment", data);
        
		kappa = getKappa("1.5");
        
        frequencies = new Frequencies();
        frequencies.initByName("data", data);
        
        instantHKY = new InstantHKY();
        instantHKY.initByName("kappa", kappa, "frequencies", frequencies);
        
        siteModel = new SiteModel();
        siteModel.initByName("gammaCategoryCount", 1, "substModel", instantHKY);
        
        pathTreeLikelihood = new PathTreeLikelihood();
        pathTreeLikelihood.initByName("PathTree", pathTree, "siteModel", siteModel);
        
        likelihood = (Distribution) pathTreeLikelihood;
        
        
        //log = new ArrayList<Plugin>();
        //log.add((Plugin)likelihood);
        //log.add((Plugin) kappa);
        //logger = new ArrayList<Logger>();
        //tmpLogger = new Logger();
        //tmpLogger.initByName("log", log);
        
        
        //logger.add(tmpLogger);
        
        // loggers
        tmpLogger = new Logger();
        tmpLogger.initByName("fileName", "testMCMC.log", "logEvery", 1, "model", likelihood, "log", likelihood, "log", kappa);
        
        // duplicate everything about the model so that it is independent of the model used in the pathTreeLikelihood calculation
		
		proposalKappa = getKappa("1.5");
        
        proposalFrequencies = new Frequencies();
        proposalFrequencies.initByName("data", data);
        
        proposalInstantHKY = new InstantHKY();
        proposalInstantHKY.initByName("kappa", proposalKappa, "frequencies", proposalFrequencies);
        
        proposalSiteModel = new SiteModel();
        proposalSiteModel.initByName("gammaCategoryCount", 1, "substModel", instantHKY);
        
        
        // operators
        allSitesPathSamplingOperator = new AllSitesPathSamplingOperator();
        allSitesPathSamplingOperator.initByName("weight", 1.0, "pathtree", pathTree, "siteModel", proposalSiteModel);

    }

	@Test
	public void testMCMCrun() throws Exception {
		
		// create  MCMC
		MCMC mcmc = new MCMC();
		//mcmc.initByName("chainLength", 100, "distribution", likelihood, "logger", logger, "operator", operatorsInput);
		
		// rightnow, use empty logger for debugging purpose
		mcmc.initByName("chainLength", 2, "distribution", likelihood, "logger", tmpLogger, "operator", allSitesPathSamplingOperator);
		// run MCMC
		mcmc.run();
	}

	
}
