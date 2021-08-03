package com.sanger.extractChains;

public class svInfo implements Comparable<svInfo>{
	private String chromLeft = null;
	private String chromRight = null;
	private Integer startLeft = 0;
	private Integer endLeft = 0;
	private Integer startRight = 0;
	private Integer endRight = 0;
	private boolean firstPos = false;
	private boolean secondPos = false;
	private String type = null;
	private Integer size = 0;
	private Boolean partChain = false;
	public svInfo(String chromLeft, Integer startLeft, Integer endLeft, String chromRight, Integer startRight, Integer endRight, boolean firstPos, boolean secondPos, String type, Integer size) {
		this.chromLeft = chromLeft;
		this.startLeft = startLeft;
		this.endLeft = endLeft;
		this.chromRight = chromRight;
		this.startRight = startRight;
		this.endRight = endRight;
		this.firstPos = firstPos;
		this.secondPos = secondPos;
		this.type = type;
		this.size = size;
	}
	public String getChromLeft() {
		return chromLeft;
	}
	public Integer getStartLeft() {
		return startLeft;
	}
	public Integer getEndLeft() {
		return endLeft;
	}
	public String getChromRight() {
		return chromRight;
	}
	public Integer getStartRight() {
		return startRight;
	}
	public Integer getEndRight() {
		return endRight;
	}
	public Boolean getFirstPos() {
		return firstPos;
	}
	public Boolean getSecondPos() {
		return secondPos;
	}
	public String getType() {
		return type;
	}
	public Integer getSize() {
		return size;
	}
	public void setChain() {
		partChain = true;
	}
	public Boolean getChain() {
		return partChain;
	}
	public void switchVariables() {
		var chrom = chromLeft;
		var start = startLeft;
		var end = endLeft;
		var first = firstPos;
		chromLeft = chromRight;
		startLeft = startRight;
		endLeft = endRight;
		firstPos = secondPos;
		chromRight = chrom;
		startRight = start;
		endRight = end;
		secondPos = first;
	}
	@Override
    public int compareTo(svInfo other) {
       return startLeft - other.getStartLeft();
    }
}