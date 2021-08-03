package com.sanger.extractChains;

import java.io.*;
import java.util.*;
import com.beust.jcommander.*;
import com.beust.jcommander.validators.PositiveInteger;
import com.sanger.extractChains.FileConverter;
import com.sanger.extractChains.FileValidator;

public class extractChains {
	
	private static String versionNumber = "0.1";
	@Parameter
	private List<String> parameters = new ArrayList<String>();
	
	@Parameter(names = "--input-bam-file", description = "Input BAM file.", required = true, converter = FileConverter.class, validateWith = FileValidator.class, order=0)
	public File input_bam_file = null;
	
	@Parameter(names = "--input-brass-bam-file", description = "Input BRASS BAM file.", required = true, converter = FileConverter.class, validateWith = FileValidator.class, order=1)
	public File input_brass_bam_file = null;
	
	@Parameter(names = "--brass-bed-file", description = "BRASS BED output file to filter.", required = true, converter = FileConverter.class, validateWith=FileValidator.class, order=2)
	public File brass_bed_file = null;
	
	@Parameter(names = "--output-file", description = "Output file to store results.", required = true, order=3)
	public String output_file = null;
	
	@Parameter(names = "--reference", description = "Reference file.", required = true, converter = FileConverter.class, validateWith=FileValidator.class, order=4)
	public File reference = null;
	
	@Parameter(names = "--width-mini-extract", description = "Mini window to extract reads (mutation_position +- width).", validateWith = PositiveInteger.class, order=5)
	public Integer extract_width_mini = 20;
	
	@Parameter(names = "--width-maxi-extract", description = "Maxi window to extract reads (mutation_position +- width).", validateWith = PositiveInteger.class, order=6)
	public Integer extract_width_maxi = 1000;
	
	@Parameter(names = "--width-minor-extract", description = "Minor window to extract reads (mutation_position +- width).", validateWith = PositiveInteger.class, order=7)
	public Integer extract_width_minor = 2000;
	
	@Parameter(names = "--width-major-extract", description = "Major window to extract reads (mutation_position +- width).", validateWith = PositiveInteger.class, order=8)
	public Integer extract_width_major = 20000;
	
	@Parameter(names = "--decoys", description = "Comma-separated list of decoy names (e.g., hs37d5).", order=9)
	public String decoys = null;
	
	@Parameter(names = {"--help","-help"}, help = true, description = "Get usage information", order=10)
	private boolean help;
	
	@Parameter(names = {"--version","-version"}, description = "Get current version", order=11)
	private Boolean version = null;
	
	public static void main(String[] args) {
		extractChains ec  = new extractChains();
		JCommander jCommander = new JCommander(ec);
		jCommander.setProgramName("extractChains");
		JCommander.newBuilder().addObject(ec).build().parse(args);
		if(ec.version != null && ec.version) {
			System.out.println("Filter BRASS output for LCM experiments: " + versionNumber);
			System.exit(0);
		}
		else if(ec.help) {
			jCommander.usage();
			System.exit(0);
		} else {
			new extractChainsCore(ec.input_bam_file, ec.input_brass_bam_file, ec.brass_bed_file, ec.output_file, ec.reference, ec.extract_width_mini, ec.extract_width_maxi, ec.extract_width_minor, ec.extract_width_major, (ec.decoys == null) ? null : new HashSet<String>(Arrays.asList(ec.decoys.split(","))));
		}
	}
}
