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
package ch.unil.genescore.vegas;

/** 
 * Interface for classes that can evaluate gene scores
 */
public abstract class GeneScoreEvaluator {
	
	
	/** Compute the gene score / p-value, return true if success */
	public abstract boolean computeScore();
	/** Get the score(s) */
	public abstract double[] getScore();

	/** Get a string representation of the results */
	public abstract String getResultsAsString();
	/** Get a header line corresponding to getResultsAsString() */
	public abstract String getResultsAsStringHeader();
	public abstract void setDataFromGeneData(GeneData Dat);		
	
	
	/** Get output to be printed on console after computing score for a gene */
	public String getConsoleOutput() { return ""; }
	
	/** Get output to be printed for genes where no score could be computed (error) */
	public String getNoScoreOutput() { return ""; }
	/** Get a header line corresponding to getNoScoreOutput() */
	public String getNoScoreOutputHeader() { return ""; }

	/** returns string defining evaluator type */
	public String getTypeString() { return ""; }

}
