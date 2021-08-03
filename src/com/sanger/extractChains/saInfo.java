package com.sanger.extractChains;

import htsjdk.samtools.*;

public class saInfo implements Comparable<saInfo>{
	public int length;
	public String chrom = null;
	public int pos = 0;
	public boolean left = false;
	public Cigar saCigar = null;
	public saInfo(int length, String chrom, int pos, boolean left, Cigar saCigar) {
		this.length = length;
		this.chrom = chrom;
		this.pos = pos;
		this.left = left;
		this.saCigar = saCigar;
	}
	
	@Override
    public int compareTo(saInfo other) {
        return other.length - this.length;
    }
	
	public int getLength() {
		return saCigar.getReferenceLength();
	}
}
