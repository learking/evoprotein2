/**
 * 
 */
package test.beast.core;


import java.util.ArrayList;

import org.junit.Before;
import org.junit.Test;

import evoprotein.proteinstructure.InputStructure;

import test.beast.evoprotein2TestCase;

import beast.core.Distribution;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Plugin;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.PathLikelihood;
import beast.evolution.likelihood.PathTreeLikelihood;
import beast.evolution.operators.PathSamplingOperator;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.InstantHKY;
import beast.evolution.substitutionmodel.ProteinCodingDNASubstModel;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.Tree;

import junit.framework.TestCase;

/**
 * @author kuangyu
 *
 */
public class MCMCourModelTest extends evoprotein2TestCase {
	
	Alignment data;
	
	Tree tree;
	PathTree pathTree;
	
	InputStructure inputStructure;
	ProteinCodingDNASubstModel ourModel;
	
	InstantHKY proposalInstantHKY;
	
	RealParameter kappa;
	RealParameter proposalKappa;
	
	Frequencies frequencies;
	Frequencies proposalFrequencies;
	
	SiteModel proposalSiteModel;
	
	PathLikelihood pathLikelihood;
	Distribution likelihood;

	Logger tmpLogger;
	ArrayList<Plugin> log;
	
	PathSamplingOperator pathSamplingOperator;

	@Before
	protected void setUp() throws Exception {
		super.setUp();
		
        data = getAlignment();
        
        tree = getTree(data);
		
        pathTree = new PathTree();
		pathTree.initByName("initial", tree, "alignment", data);
        
		kappa = getKappa("1.5");
        
        frequencies = new Frequencies();
        frequencies.initByName("data", data);
        
        //our model
        inputStructure = new InputStructure();
        // inputStructure.initByName();
        ourModel = new ProteinCodingDNASubstModel();
        // ourModel.initByName("kappa", kappa, "frequencies", frequencies, "inputStructure", );
        
        pathLikelihood = new PathLikelihood();
        //pathLikelihood.initByName("PathTree", pathTree, "ourModel", ourModel);
        
        likelihood = (Distribution) pathLikelihood;
        
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
        proposalSiteModel.initByName("gammaCategoryCount", 1, "substModel", proposalInstantHKY);
        
        
        // operators
        pathSamplingOperator = new PathSamplingOperator();
        pathSamplingOperator.initByName("weight", 1.0, "pathtree", pathTree, "siteModel", proposalSiteModel);
		
	}

	@Test
	public void testMCMCourModelRun() throws Exception {
		// create  MCMC
		MCMC mcmc = new MCMC();
		// right now, use empty logger for debugging purpose
		mcmc.initByName("chainLength", 1000, "distribution", likelihood, "logger", tmpLogger, "operator", pathSamplingOperator);
		// run MCMC
		mcmc.run();
	}

}
