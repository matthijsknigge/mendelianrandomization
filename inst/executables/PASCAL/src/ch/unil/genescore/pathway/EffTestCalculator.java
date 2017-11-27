/*******************************************************************************
 * Copyright (c) 2015 David Lamparter, Daniel Marbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *******************************************************************************/
package ch.unil.genescore.pathway;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import ch.unil.genescore.gene.Gene;

/**
 * Calculates the effective number of Tests for the given pathway library.
 * Atm only implements chi2-statistic
 */
public class EffTestCalculator {
	/**genes_: genes that contain random scores. these scores are sampled.
	 * */
	private HashMap<String,Gene> genes_ = null;
	private HashMap<String,Integer> indexMap_ = new  HashMap<String,Integer>();
			//null;
	private double[] geneVals_ = null; 
	/**geneSets_: geneSets_ that are used to get ids:
	 * these are the genes that have been calculated: Do not change internally.
	 * */	
	private ArrayList<GeneSet> geneSets_ = null;
	private static ChiSquaredDistribution chiSquared1df_ = new ChiSquaredDistribution(1);
	
	
	public void addGenes(Collection<? extends Gene> genes){
		for (Gene gene : genes){
			Gene fake = new Gene(gene.id_);
			fake.setChi2Stat(chiSquared1df_.sample());
			genes_.put(gene.id_,fake);			
		}	
	}
	public void addMetaGenes(HashMap<String, MetaGene> metaGenes){		
		for (Entry<String, MetaGene> e : metaGenes.entrySet()){
			Gene g =  e.getValue();
			Gene fake = new Gene(g.id_);
			fake.setChi2Stat(chiSquared1df_.sample());
			genes_.put(g.id_,fake);			
		}
	}
	
	public void initializeGeneVals(){
		geneVals_ = new double[indexMap_.size()];
	}
	
	
	public void addGeneIds(Collection<? extends Gene> genes){
		for (Gene gene : genes){
			indexMap_.put(gene.id_,indexMap_.size());						
		}	
	}
	public void addMetaGeneIds(HashMap<String, MetaGene> metaGenes){		
		for (Entry<String, MetaGene> e : metaGenes.entrySet()){
			Gene gene =  e.getValue();
			indexMap_.put(gene.id_,indexMap_.size());	
		}
	}
	protected void calcGeneVals(){
		
		for (int i=0 ; i< geneVals_.length; i++){
			geneVals_[i]=chiSquared1df_.sample();				
		}
	}
	
	

	
	private double getMinScores(double[] score){
		double min=1;
		for (double d : score){
			if (d<min){
				min=d;
			}
		}
		return(min);
	}
	
	public double calcMinVal(){
		calcGeneVals();
		double min = getMinScores(calcSetScores());
		return min;
	}
	
	private int[] getIndexSet(ArrayList<String> ids){
		int[] indSet= new int[ids.size()];		
		for (int i=0 ; i< ids.size();i++){
			String id = ids.get(i);
			if(!indexMap_.containsKey(id)){
				throw new RuntimeException("Not all ids present");
			}		
			indSet[i]=indexMap_.get(id);		
		}
		return(indSet);
	}
	
	double getScoreSum(int[] indSet){
		double scoreSum=0;
		for (int i=0; i<indSet.length; i++){
		scoreSum = scoreSum + geneVals_[indSet[i]];		
		}
		return(scoreSum);
	}
	
	private double[] calcSetScores(){
		ChiSquaredDistribution chiSquaredNdf = null;
		double[] setScores = new double[geneSets_.size()];
		for (int i=0 ; i<geneSets_.size() ;i++){
			ArrayList<String> gIds = geneSets_.get(i).getGeneIds();
			int[] indSet = getIndexSet(gIds);
			double scoreSum = getScoreSum(indSet);
			chiSquaredNdf = new ChiSquaredDistribution(indSet.length);
			setScores[i] = 1-chiSquaredNdf.cumulativeProbability(scoreSum); 	
		}		
		return(setScores);
	}
			
	public void recalcFakeScores(){
		Set<Entry<String,Gene>> gg = genes_.entrySet();
		for (Entry<String,Gene> g: gg){
			g.getValue().setChi2Stat(chiSquared1df_.sample());
		}		
	}
			
	public void setGeneSets(ArrayList<GeneSet> geneSets){
		geneSets_= geneSets;		
	}
	
	
	
	
}
