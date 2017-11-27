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
package ch.unil.genescore.vegas.test;

import java.util.ArrayList;
import java.util.Collections;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.vegas.OverlappedGenomicElement;
import static org.junit.Assert.*;


import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import ch.unil.genescore.vegas.OverlappedCollectionStream;



/**
 * Unit tests for GenewiseSnpWeightsTest
 */
public class OverlappedCollectionStreamTest {

	
	static double delta;
	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {	
		delta=1E-14;
	}

	@AfterClass
	public static void testCleanup() { }
	
	@Rule
    public ExpectedException thrown= ExpectedException.none();
	
	// ============================================================================
	// TESTS

	/** Test PathwayMain.run() */
	@Test
	public void  RunCollectionStream1Test() {
	// test standard application
	String chr1="chr19";
	String chr2="chr22";
	Gene g1 = new Gene("g1");
	Gene g2 = new Gene("g2");
	Gene g3 = new Gene("g3");
	g1.setPosition(chr1, 1000, 2000, true);
	g2.setPosition(chr1, 1900, 3000, true);
	g3.setPosition(chr2, 1000, 2000, true);
	
	ArrayList<OverlappedGenomicElement> gArray = new ArrayList<OverlappedGenomicElement>(3);
	gArray.add(new OverlappedGenomicElement(g1));
	gArray.add(new OverlappedGenomicElement(g2));
	gArray.add(new OverlappedGenomicElement(g3));
	Collections.sort(gArray);
	OverlappedCollectionStream str = new OverlappedCollectionStream(gArray);
	OverlappedGenomicElement el1 = str.getNextAsOverlappedGenomicElement();
	OverlappedGenomicElement el2 = str.getNextAsOverlappedGenomicElement();
	OverlappedGenomicElement el3 = str.getNextAsOverlappedGenomicElement();
	
	assertTrue(el1.getMainElement()==g1);
	assertTrue(el2.getMainElement()==g2);
	assertTrue(el3.getMainElement()==g3);	
	}	
	@Test
	public void  RunCollectionStreamTest2() {
		//Test closed stream exception
	String chr1="chr19";
	Gene g1 = new Gene("g1");
	g1.setPosition(chr1, 1000, 2000, true);	
	ArrayList<OverlappedGenomicElement> gArray = new ArrayList<OverlappedGenomicElement>(3);
	gArray.add(new OverlappedGenomicElement(g1));	
	Collections.sort(gArray);
	OverlappedCollectionStream str = new OverlappedCollectionStream(gArray);
	OverlappedGenomicElement el1 = str.getNextAsOverlappedGenomicElement();
	thrown.expect(RuntimeException.class);
	OverlappedGenomicElement el2 = str.getNextAsOverlappedGenomicElement();
	}	
	@Test
	public void  RunCollectionStream3Test() {
		// test unsorted exception
		String chr1="chr2";
		String chr2="chr12";
		Gene g1 = new Gene("g1");
		Gene g2 = new Gene("g2");
		Gene g3 = new Gene("g3");
		g1.setPosition(chr1, 1000, 2000, true);
		g2.setPosition(chr1, 1900, 3000, true);
		g3.setPosition(chr2, 1000, 2000, true);
		
		ArrayList<OverlappedGenomicElement> gArray = new ArrayList<OverlappedGenomicElement>(3);
		
		gArray.add(new OverlappedGenomicElement(g1));
		gArray.add(new OverlappedGenomicElement(g2));
		gArray.add(new OverlappedGenomicElement(g3));
		
		OverlappedCollectionStream str = new OverlappedCollectionStream(gArray);
		thrown.expect(RuntimeException.class);
		OverlappedGenomicElement el1 = str.getNextAsOverlappedGenomicElement();
		OverlappedGenomicElement el2 = str.getNextAsOverlappedGenomicElement();
		OverlappedGenomicElement el3 = str.getNextAsOverlappedGenomicElement();			
		}	
}