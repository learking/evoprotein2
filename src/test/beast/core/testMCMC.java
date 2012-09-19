/**
 * 
 */
package test.beast.core;


import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.Plugin;
import beast.core.parameter.Parameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.PathTreeLikelihood;
import beast.evolution.operators.PathSamplingOperator;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.InstantHKY;
import beast.evolution.tree.PathTree;
import beast.evolution.tree.Tree;

import test.beast.evoprotein2TestCase;

import junit.framework.TestCase;

/**
 * @author kuangyu
 *
 */
public class testMCMC extends evoprotein2TestCase {

	Alignment data;
	Tree tree;
	PathTree pathTree;
	InstantHKY instantHKY;
	RealParameter kappa;
	Frequencies frequences;
	SiteModel siteModel;
	PathTreeLikelihood pathTreeLikelihood;
	Distribution likelihood;

	ArrayList<Logger> logger;
	Logger tmpLogger;
	ArrayList<Plugin> log;
	ArrayList<Operator> operatorsInput;
	Operator tmpOperator;
	
	PathSamplingOperator pathSamplingOperator;
	
    @Override
    protected void setUp() throws Exception {
        super.setUp();
        
        data = getAlignment();
        
        tree = getTree(data);
		
        pathTree = new PathTree();
		pathTree.initByName("initial", tree, "alignment", data);
        
		kappa = getKappa();
        
        frequences = new Frequencies();
        frequences.initByName("data", data);
        
        instantHKY = new InstantHKY();
        instantHKY.initByName("kappa", kappa, "frequencies", frequences);
        
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
        /*
        List<Logger> m_logger = new ArrayList<Logger>();
        tmpLogger = new Logger();
        tmpLogger.initByName("log", likelihood);
        m_logger.add(tmpLogger);
        tmpLogger = new Logger();
        tmpLogger.initByName("log", kappa);
        m_logger.add(tmpLogger);
        */
        // why type mismatch?
        
        // operators
        operatorsInput = new ArrayList<Operator>();
        pathSamplingOperator = new PathSamplingOperator();
        pathSamplingOperator.initByName("weight", 1.0, "pathtree", pathTree, "siteModel", siteModel);
        operatorsInput.add(pathSamplingOperator);
        // why type mismatch?
    }

	@Test
	public void test() throws Exception {
		
		// create  MCMC
		MCMC mcmc = new MCMC();
		//mcmc.initByName("chainLength", 100, "distribution", likelihood, "logger", logger, "operator", operatorsInput);
		
		// rightnow, use empty logger for debugging purpose
		mcmc.initByName("chainLength", 100, "distribution", likelihood, "logger", logger, "operator", pathSamplingOperator);
		// run MCMC
		mcmc.run();
		
		System.out.println("Setup is a success!");
	}

	
}
