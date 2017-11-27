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
package ch.unil.genescore.gene;

import java.io.File;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Map.Entry;

import ch.unil.genescore.main.Pascal;
import ch.unil.gpsutils.FileExport;
import ch.unil.gpsutils.FileParser;


/**
 * Gene annotation, provides functionality to load IDs and coordinates of genes
 */
abstract public class GeneAnnotation {

	/** The file with the genome annotation */
	protected File annotationFile_ = null;
	
	/** User specified set of genes to be loaded (leave empty for all genes, boolean indicates if gene was found in the annotation) */
	protected HashMap<String, Boolean> genesToBeLoaded_ = null;
	/** Load only genes from this chromosome (leave empty for all chromosomes) */
	protected String chromosomeToBeLoaded_ = null;
	/** Ignore X and Y chromosomes */
	protected boolean ignoreAllosomes_ = false;
	/** Ignore mitochondiral genes */
	protected boolean ignoreChrM_ = false;
	/** Flag indicates whether only protein coding genes should be loaded */
	protected boolean loadOnlyProteinCoding_ = false;
	
	/** The genes that were loaded from the annotation file */
	protected LinkedHashMap<String, Gene> genes_ = null;
	/** Set of genes that were excluded */
	protected HashSet<String> excludedGenes_ = null;
	
	
	// ============================================================================
	// STATIC METHODS

	/** Create an instance of the appropriate subclass as specified in Settings.genomeAnnotation_ */
	public static GeneAnnotation createAnnotationInstance() {
		
		// TODO allow annotation to be loaded from a bed file, create appropriate subclass
		if (Pascal.set.genomeAnnotation_.equals("gencode"))
			return new GeneAnnotationGencode();
		else if (Pascal.set.genomeAnnotation_.equals("ucsc"))
			return new GeneAnnotationUcsc();
		else if (Pascal.set.genomeAnnotation_.equals("bed"))
			return new GeneAnnotationBed();
			throw new RuntimeException("Settings.genomeAnnotation_ must be 'gencode','ucsc' or 'bed'.");
	}
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public GeneAnnotation(File annotationFile) {
		
		annotationFile_ = annotationFile; 
		chromosomeToBeLoaded_ = Pascal.set.chromosome_;
		ignoreAllosomes_ = Pascal.set.ignoreAllosomes_;
		ignoreChrM_ = Pascal.set.ignoreChrM_;
		loadOnlyProteinCoding_ = Pascal.set.loadOnlyProteinCodingGenes_;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Load the specified genes (genesToBeLoaded_) from the annotation file */
	public LinkedList<Gene> loadAnnotation(File genesToBeLoadedFile) {
		
		// Load the set of genes to be considered
		loadGenesToBeLoaded(genesToBeLoadedFile);
		// Load the specified genes from the annotation
		loadAnnotation();
		
		// If not all specified genes were found, print warning
		if (genesToBeLoaded_ != null && genesToBeLoaded_.size() != genes_.size()) {
			int numNotFound = genesToBeLoaded_.size() - genes_.size();
			
			String genesNotFound = "";
			for (Entry<String,Boolean> entry : genesToBeLoaded_.entrySet())
				if (!entry.getValue())
					genesNotFound += entry.getKey() + " ";

			Pascal.println("   - " + numNotFound + " genes were not found in the annotation: " + genesNotFound);
		}
		//Don't need the hashmap put into list
		LinkedList<Gene> genes = new LinkedList<Gene>();
		for(Gene gene : genes_.values()){
			genes.add(gene);
		}
		Collections.sort(genes);
		return genes;
	}

	// ----------------------------------------------------------------------------

	/** Load the specified genes (genesToBeLoaded_) from the annotation file */
	abstract public LinkedHashMap<String, Gene> loadAnnotation();

	
	// ----------------------------------------------------------------------------

	/** Load the specified set of genes (genesToBeLoaded_) */
	public void loadGenesToBeLoaded(File file) {
				
		// Return if no gene file was specified (all genes from the annotation will be loaded)
		if (file == null)
			return;

		genesToBeLoaded_ = new HashMap<String, Boolean>();
		FileParser parser = new FileParser(Pascal.log, file);
		
		while (true) {
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;
			if (nextLine.length != 1)
				parser.error("Expected one column (gene ID)");
			
			genesToBeLoaded_.put(nextLine[0].toUpperCase(), false);
		}
		parser.close();
	}


	// ----------------------------------------------------------------------------

	/** Write the genes to a file with gene id, symbol and position */
	public void writeGeneList(File file) {
		
		FileExport writer = new FileExport(Pascal.log, file);
		String prevChr = null;
		int prevStart = -1;
		
		String header = "geneId\tgeneSymbol\tchromosome\tstart\tend\tstrand";
		writer.println(header);
		
		for (Gene gene : genes_.values()) {
			if (prevStart == -1 || prevChr != gene.chr_) {
				prevStart = gene.start_;
				prevChr = gene.chr_;
			}
			
			if (gene.start_ < prevStart)
				Pascal.error("Genes are not ordered by genomic position");
			
			String nextLine = gene.id_ + "\t" + gene.symbol_ + "\t" + 
					gene.chr_ + "\t" + gene.start_ + "\t" + gene.end_ + "\t" + 
					(gene.posStrand_ ? "+" : "-");
			
			writer.println(nextLine);
		}
		writer.close();
	}
	
	
	// ----------------------------------------------------------------------------

	/** Return true if this chromosome should be skipped */
	public boolean skipChromosome(String chr) {
		
		boolean notSpecifiedChrom = chromosomeToBeLoaded_ != null && !chromosomeToBeLoaded_.equals("") && !chr.equals(chromosomeToBeLoaded_);
		boolean ignoreChrom = (ignoreAllosomes_ && (chr.equals("chrX") || chr.equals("chrY"))) ||
				(ignoreChrM_ && chr.equals("chrM"));
		boolean isValid = Chromosome.isValidId(chr);
		
		return !isValid || notSpecifiedChrom || ignoreChrom;
	}
	

	// ----------------------------------------------------------------------------

	/** Remove the genes in the given file from the loaded set of genes */
	public int removeGenes(File excludedGenesFile) {
		
		excludedGenes_ = loadGeneList(excludedGenesFile);
		int numRemoved = 0;
		
		for (String gene : excludedGenes_) {
			Gene g = genes_.remove(gene);
			if (g != null)
				numRemoved++;
		}
		return numRemoved;
	}

	
	// ----------------------------------------------------------------------------
	
	/** Load a list of genes (one line per gene, id in first column) */
	public static HashSet<String> loadGeneList(File excludedGenesFile) {
		
		HashSet<String> excludedGenes = new HashSet<String>();
		if (excludedGenesFile == null)
			return excludedGenes;
		
		FileParser parser = new FileParser(Pascal.log, excludedGenesFile);

		// Parse header
		String[] header = parser.readLine();
		if (!header[0].equals("gene_id"))
			Pascal.error("Expected header line with first field 'gene_id' (tab-separated)");
		
		// First line
		while (true) {
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;
			
			// Gene id
			String id = nextLine[0];
			id = GeneAnnotationGencode.removeEnsemblVersion(id);
			excludedGenes.add(id);
		}
		parser.close();
		
		return excludedGenes;
	}

	
	/** Get name of genome annotation for console output ('UCSC known gene', 'GENCODE genes') */
	public static String getAnnotationName() {
		
		if (Pascal.set.genomeAnnotation_.equalsIgnoreCase("ucsc"))
			return "UCSC known genes";
		else if (Pascal.set.genomeAnnotation_.equalsIgnoreCase("gencode"))
			return "GENCODE genes";
		else
			throw new RuntimeException("Unknown genomeAnnotation:" + Pascal.set.genomeAnnotation_);
	}
	
	
	// ============================================================================
	// PROTECTED METHODS

	/** Check if the given gene ID or symbol are listed in genesToBeLoaded_, if yes return true and set the corresponding entry true */
	protected boolean checkGenesToBeLoaded(String geneId, String geneName) {
	
		if (genesToBeLoaded_ == null)
			return true;
		
		// Check if this gene is part of the specified gene set
		String specifiedGene = null;
		if (genesToBeLoaded_.containsKey(geneId))
			specifiedGene = geneId;
		else if (genesToBeLoaded_.containsKey(geneName))
			specifiedGene = geneId;

		// Flag this gene as found
		if (specifiedGene != null) {
			genesToBeLoaded_.put(specifiedGene, true);
			return true;

		} else {
			return false;
		}
	}
	
}
