package com.sanger.extractChains;

public class chainResults {
	public int fbi = 0;
	public int fbi_chain = 0;
	public int fbi_broken = 0;
	public int rte = 0;
	public int cli = 0;
	public chainResults() {}
	
	public void increaseFbi() {
		fbi++;
	}
	public void increaseFbiChain() {
		fbi_chain++;
	}
	public void increaseFbiChainBroken() {
		fbi_broken++;
	}
	public void increaseRte() {
		rte++;
	}
	public void increaseCli() {
		cli++;
	}
}
