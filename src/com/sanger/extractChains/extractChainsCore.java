package com.sanger.extractChains;

import java.io.*;
import java.util.*;

import com.sanger.intervalTree.Interval1D;
import com.sanger.intervalTree.IntervalST;

import htsjdk.samtools.*;
import htsjdk.samtools.cram.build.CramIO;


public class extractChainsCore {
	private int nano = 5;
	private ArrayList<svInfo> svList = null;
	public extractChainsCore(File input_bam, File input_brass, File bed, String output, File reference, int mini, int maxi, int minor, int major, HashSet<String> decoys) {
		System.out.println("Retrieving SVs...");
		var buffer = new StringBuffer();
		svList = getStructuralVariants(bed);
		var bpTree = new IntervalST<bpInfo>();
		var results = firstGetEvents(svList, input_bam, mini, maxi, minor, major, decoys, new chainResults(), reference, false, bpTree, buffer);
		svList = reorgList(svList);
		Collections.sort(svList);
		results = firstGetEvents(svList, input_bam, mini, maxi, minor, major, decoys, results, reference, true, bpTree, buffer);
		writeResults(output, buffer);
		System.out.printf("FBI: %d - FBI broken: %d - FBI full chain: %d - Retrotransposon insertion: %d - Classical inversion: %d%n", results.fbi, results.fbi_broken, results.fbi_chain, results.rte, results.cli);
	}
	private void writeResults(String output, StringBuffer buffer) {
		try {
			var out = new BufferedWriter(new FileWriter(output));
			out.write(String.join("\t", "First chrom left", "First start left", "First chrom right", "First start right", "Sec chrom left", "Sec start left", "Sec chrom right", "Sec start right", "Class", "Chainlength 1", "Chainlength 2") + "\n" + buffer.toString());
			out.flush();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-3);
		}
		
	}
	private ArrayList<svInfo> reorgList(ArrayList<svInfo> svList) {
		svInfo current = null;
		for(int i = 0; i < svList.size(); i++) {
			current = svList.get(i);
			current.switchVariables();
			svList.set(i, current);
		}
		return svList;
	}
	private chainResults firstGetEvents(ArrayList<svInfo> svList, File bam, int mini, int maxi, int minor, int major, HashSet<String> decoys, chainResults results, File reference, boolean sorted, IntervalST<bpInfo> bpTree, StringBuffer buffer) {
		int i,j;
		svInfo current = null;
		svInfo next = null;
		svloop:
		for(i = 0; i < svList.size(); i++) {
			var range = new ArrayList<svInfo>(100); 
			current = svList.get(i);
			if(current.getChain())
				continue svloop;
			if(partOfTree(bpTree, current.getChromLeft(), current.getStartLeft(), current.getEndLeft())) {
				//System.out.printf("%s - %d - %d - Analysed before%n", current.getChromLeft(), current.getStartLeft(), current.getEndLeft());
				continue svloop;
			}
			for(j = i+1; j < svList.size(); j++) {
				next = svList.get(j);
				if(next.getChain())
					continue;
				if(current.getChromLeft().equals(next.getChromLeft()) && (next.getStartLeft() - current.getStartLeft()) <= major)
					range.add(next);
				else if(sorted && (next.getStartLeft() - current.getStartLeft()) <= major) {}
				else
					break;
			}
			if(current.getType().equals("inversion") && current.getFirstPos() != current.getSecondPos()) {
				for(j = 0; j < range.size(); j++) {
					next = range.get(j);
					if(next.getType().equals("inversion") && Math.abs(current.getStartLeft() - next.getStartLeft()) <= maxi && Math.abs(current.getStartRight() - next.getStartRight()) <= maxi && current.getFirstPos() != next.getFirstPos()) {
						current.setChain();
						next.setChain();
						results.increaseCli();
						buffer.append(String.join("\t", current.getChromLeft(), current.getStartLeft().toString(), current.getChromRight(), current.getStartRight().toString(), next.getChromLeft(), next.getStartLeft().toString(), next.getChromRight(), next.getStartRight().toString(), "Classical inversion", "1", "1") + "\n");
						continue svloop;
					}
				}
				if(current.getSize() <= minor) {
					current.setChain();
					buffer.append(String.join("\t", current.getChromLeft(), current.getStartLeft().toString(), current.getChromRight(), current.getStartRight().toString(), "-", "-", "-", "-", "FBI single event", "1", "1") + "\n");
					results.increaseFbi();
				}
			} else if(current.getType().equals("translocation")) {
				String len = null;
				for(j = 0; j < range.size(); j++) {
					next = range.get(j);
					if(current.getType().equals(next.getType()) && rteEvent(current, next, bam, reference, mini)) {
						current.setChain();
						next.setChain();
						results.increaseRte();
						buffer.append(String.join("\t", current.getChromLeft(), current.getStartLeft().toString(), current.getChromRight(), current.getStartRight().toString(), next.getChromLeft(), next.getStartLeft().toString(), next.getChromRight(), next.getStartRight().toString(), "RTE event", "1", "1") + "\n");
						continue svloop;
					}
					else if(current.getType().equals(next.getType()) && current.getFirstPos() == next.getFirstPos()) {
						var firstChain = chaining(current, next, bam, reference, mini, maxi, decoys, sorted, bpTree);
						System.out.println("Stop 1");
						if(firstChain.hindsightRTE) {
							current.setChain();
							next.setChain();
							results.increaseRte();
							System.out.println("Hindsight RTE");
							buffer.append(String.join("\t", current.getChromLeft(), current.getStartLeft().toString(), current.getChromRight(), current.getStartRight().toString(), next.getChromLeft(), next.getStartLeft().toString(), next.getChromRight(), next.getStartRight().toString(), "RTE event", "1", "1") + "\n");
							continue svloop;
						}
						else if(!firstChain.broken && firstChain.chains > 1) {
							current.setChain();
							next.setChain();
							results.increaseFbiChain();
							System.out.println("Chain completed");
							buffer.append(String.join("\t", current.getChromLeft(), current.getStartLeft().toString(), current.getChromRight(), current.getStartRight().toString(), next.getChromLeft(), next.getStartLeft().toString(), next.getChromRight(), next.getStartRight().toString(), "FBI chain complete", firstChain.chains.toString(), firstChain.chains.toString()) + "\n");
							continue svloop;
						} else if(firstChain.broken && firstChain.chains >= 2) {
							current.setChain();
							next.setChain();
							results.increaseFbiChainBroken();
							System.out.println("Chain broken");
							var secondChain = chaining(next, current, bam, reference, mini, maxi, decoys, sorted, bpTree);
							buffer.append(String.join("\t", current.getChromLeft(), current.getStartLeft().toString(), current.getChromRight(), current.getStartRight().toString(), next.getChromLeft(), next.getStartLeft().toString(), next.getChromRight(), next.getStartRight().toString(), "FBI chain broken", firstChain.chains.toString(), secondChain.chains.toString()) + "\n");
							continue svloop;
						} else
							len = firstChain.chains.toString();
						firstChain = chaining(next, current, bam, reference, mini, maxi, decoys, sorted, bpTree);
						System.out.println("Stop 2");
						if(!firstChain.broken && firstChain.chains > 1) {
							current.setChain();
							next.setChain();
							results.increaseFbiChain();
							System.out.println("Chain completed");
							buffer.append(String.join("\t", current.getChromLeft(), current.getStartLeft().toString(), current.getChromRight(), current.getStartRight().toString(), next.getChromLeft(), next.getStartLeft().toString(), next.getChromRight(), next.getStartRight().toString(), "FBI chain complete", firstChain.chains.toString(), firstChain.chains.toString()) + "\n");
							continue svloop;
						} else if(firstChain.chains >= 2) {
							current.setChain();
							next.setChain();
							results.increaseFbiChainBroken();
							System.out.println("Chain broken");
							buffer.append(String.join("\t", current.getChromLeft(), current.getStartLeft().toString(), current.getChromRight(), current.getStartRight().toString(), next.getChromLeft(), next.getStartLeft().toString(), next.getChromRight(), next.getStartRight().toString(), "FBI chain broken", len, firstChain.chains.toString()) + "\n");
							continue svloop;
						}
					}
				}
			}
		}
		return results;
	}
	private boolean partOfTree(IntervalST<bpInfo> bpTree, String chrom, int start, int end) {
		var bpIterator = bpTree.searchAll(new Interval1D(start - nano, end + nano)).iterator();
		bpInfo bp = null;
		while(bpIterator.hasNext()) {
			bp = bpTree.get(bpIterator.next());
			if(chrom.equals(bp.chrom))
				return true;
		}
		return false;
	}
	private void updateTree(IntervalST<bpInfo> bpTree, String chrom, int start, int end) {
		boolean found = false;
		var it = bpTree.searchAll(new Interval1D(start, end)).iterator();
		bpInfo bp = null;
		while(it.hasNext()) {
			bp = bpTree.get(it.next());
			if(bp.chrom.equals(chrom)) {
				found = true;
				break;
			}
		}
		if(!found)
			bpTree.put(new Interval1D(start, end), new bpInfo(chrom, start)); 
	}
	private chainInfo chaining(svInfo current, svInfo next, File bam, File reference, int mini, int maxi, HashSet<String> decoys, boolean sorted, IntervalST<bpInfo> bpTree) {
		System.out.println("Entered chaining analysis...");
		var cSet = new HashSet<String>(100, 0.9999f);
		var chain = new chainInfo();
		SAMRecordIterator it = null;
		SamReader inputSam = null;
		boolean go = true;
		Cigar cig = null;
		Cigar saCigar = null;
		CigarElement tel = null;
		CigarElement saLeft = null;
		CigarElement saRight = null;
		String[] tokens = null;
		String[] saItems = null;
		String[] saInfo = null;
		String sa = null;
		String chrom = current.getChromLeft();
		Integer start = current.getStartLeft();
		Integer end = current.getEndLeft();
		int pos = -1;
		if(bam.toString().endsWith(".cram") && CramIO.checkHeaderAndEOF(bam))
			inputSam = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).referenceSequence(reference).open(SamInputResource.of(bam));
		else
			inputSam = SamReaderFactory.make().enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX).validationStringency(ValidationStringency.LENIENT).samRecordFactory(DefaultSAMRecordFactory.getInstance()).open(bam);
		while(go) {
			updateTree(bpTree, chrom, start, end);
			if(cSet.contains(chrom + "_" + start.toString())) {
				System.out.println("Ouroboros achieved");
				chain.chains = 0;
				return chain;
			} else
				cSet.add(chrom + "_" + start.toString());
			System.out.printf("%s - %s - %d - %d%n", (!sorted) ? "Forward" : "Reverse", chrom, start, end);
			var saList = new ArrayList<saInfo>(200);
			if(hindsightRTE(inputSam, chrom, start, end)) {
				chain.hindsightRTE = true;
				return chain;
			}
			it = inputSam.query(chrom, start - mini, end + mini, false);
			SAMRecord currentRecord = null;
			while(it.hasNext()) {
				currentRecord = it.next();
				sa = currentRecord.getStringAttribute("SA");
				if(sa == null || sa.equals(""))
					continue;
				if(sa.split(";").length == 2) {
					saItems = sa.split(";");
					int point = 0;
					int diff = 0;
					for(int j = 0; j < saItems.length; j++) {
						saInfo = saItems[j].split(",");
						pos = Integer.parseInt(saInfo[1]);
						if(!(saInfo[0].equals(chrom) && (Math.abs(start - pos) <= mini || Math.abs(start - (pos + TextCigarCodec.decode(saInfo[3]).getReferenceLength())) <= mini))) {
							point = j;
							diff++;
						}
					}
					if(diff == 1)
						sa = saItems[point];
					else {
						boolean found = false;
						for(int j = 0; j < saItems.length; j++) {
							saInfo = saItems[j].split(",");
							pos = Integer.parseInt(saInfo[1]);
							saCigar = TextCigarCodec.decode(saInfo[3]);
							if(saInfo[0].equals(next.getChromLeft()) && (Math.abs(pos - next.getStartLeft()) <= maxi || Math.abs((pos + saCigar.getReferenceLength()) - next.getStartLeft()) <= maxi)) {
								sa = saItems[j];
								found = true;
								break;
							}
						}
						if(!found)
							continue;
					}
				} else if(sa.split(";").length > 2) {
					//System.out.println("We are in trouble");
					continue;
				}
				tokens = sa.split(",");
				cig = currentRecord.getCigar();
				saCigar = TextCigarCodec.decode(tokens[3]);
				saLeft = saCigar.getFirstCigarElement();
				saRight = saCigar.getLastCigarElement();
				if(Math.abs(currentRecord.getAlignmentStart() - start) <= nano) {
					tel = cig.getFirstCigarElement();
					if(tel.getOperator().toString().equals("S"))
						saList.add(new saInfo(tel.getLength(), tokens[0], Integer.parseInt(tokens[1]), (saCigar.isLeftClipped() && saCigar.isRightClipped()) ? (saLeft.getLength() >= saRight.getLength()) : saCigar.isLeftClipped(), saCigar));
				} else if(Math.abs(currentRecord.getAlignmentEnd() - start) <= nano) {
					tel = cig.getLastCigarElement();
					if(tel.getOperator().toString().equals("S"))
						saList.add(new saInfo(tel.getLength(), tokens[0], Integer.parseInt(tokens[1]), (saCigar.isLeftClipped() && saCigar.isRightClipped()) ? (saLeft.getLength() >= saRight.getLength()) : saCigar.isLeftClipped(), saCigar));
				}
			}
			it.close();
			if(saList.size() < 5) {
				chain.broken = true;
				break;
			}
			Collections.sort(saList);
			int count = 0;
			saInfo first = saList.get(0);
			for(int i = 1; i < 5; i++) {
				if(first.chrom.equals(saList.get(i).chrom) && first.saCigar.isLeftClipped() && Math.abs(first.pos - saList.get(i).pos) <= nano)
					count++;
				else if(first.chrom.equals(saList.get(i).chrom) && first.saCigar.isRightClipped() && Math.abs((first.pos + first.saCigar.getReferenceLength()) - (saList.get(i).pos + saList.get(i).saCigar.getReferenceLength())) <= nano)
					count++;
			}
			if(count <= 3 || decoys.contains(first.chrom)) {
				chain.broken = true;
				break;
			} else
				chain.chains = chain.chains + 1;
			if(first.chrom.equals(next.getChromLeft()) && (Math.abs(first.pos - next.getStartLeft()) <= maxi || Math.abs((first.pos + first.saCigar.getReferenceLength()) - next.getStartLeft()) <= maxi))
				return chain;
			updateTree(bpTree, first.chrom, (first.left) ? first.pos : first.pos + first.getLength(), (first.left) ? first.pos + 1 : first.pos + first.getLength() + 1);
			it = (first.left) ? inputSam.query(first.chrom, first.pos, first.pos + maxi, false) : inputSam.query(first.chrom, first.pos + first.getLength() - maxi, first.pos + first.getLength(), false);
			pos = (first.left) ? first.pos : first.pos + first.getLength() - maxi;
			var hist = new int[maxi];
			int firstElement = 0;
			int runningSum = 0;
			int max = 0;
			int winner = -1;
			while(it.hasNext()) {
				currentRecord = it.next();
				if(first.left && currentRecord.getCigar().isRightClipped() && (currentRecord.getAlignmentEnd() - pos) < maxi) {
					//System.out.printf("Right: %d%n", currentRecord.getAlignmentEnd() - pos);
					hist[currentRecord.getAlignmentEnd() - pos] += 1;
				} else if(!first.left && currentRecord.getCigar().isLeftClipped() && (currentRecord.getAlignmentStart() - pos) >= 0 && (currentRecord.getAlignmentStart() - pos) < maxi) {
					//System.out.printf("Left: %d%n", currentRecord.getAlignmentStart() - pos);
					hist[currentRecord.getAlignmentStart() - pos] += 1;
				}
			}
			it.close();
			for(int i = 0; i < nano; i++)
				runningSum += hist[i];
			firstElement = 0;
			for(int i = 1; i < (hist.length - nano); i++) {
				runningSum -= firstElement;
				runningSum += hist[i + nano - 1];
				firstElement = hist[i];
				if(runningSum > max) {
					max = runningSum;
					winner = i;
				}
			}
			if(max >= 5) {
				chrom = first.chrom;
				start = end = pos + winner;
			} else {
				chain.broken = true;
				return chain;
			}
			
		}
		return chain;
	}
	private boolean hindsightRTE(SamReader input, String chrom, int start, int end) {
		var it = input.query(chrom, start - nano, end + nano, false);
		double count = 0, total = 0;
		CigarElement tel = null;
		SAMRecord current = null;
		Cigar currentCigar = null;
		CigarElement currentTel = null;
		while(it.hasNext()) {
			current = it.next();
			double polyA = 0;
			double polyT = 0;
			byte[] seqRead = current.getReadBases();
			int trackerSeq = 0;
			currentCigar = current.getCigar();
			if(Math.abs(current.getAlignmentStart() - start) <= nano && currentCigar.isLeftClipped() && currentCigar.getFirstCigarElement().getLength() >= 20 && !currentCigar.getFirstCigarElement().getOperator().toString().equals("H")) {
				total++;
				tel = currentCigar.getFirstCigarElement();
				for(int i = 0; i < tel.getLength(); i++) {
					if(seqRead[trackerSeq] == 'A')
						polyA++;
					else if (seqRead[trackerSeq] == 'T')
						polyT++;
					trackerSeq++;
				}
			} else if(Math.abs(current.getAlignmentEnd() - start) <= nano && currentCigar.isRightClipped() && currentCigar.getLastCigarElement().getLength() >= 20 && !currentCigar.getLastCigarElement().getOperator().toString().equals("H")) {
				total++;
				tel = currentCigar.getLastCigarElement();
				var elList = currentCigar.getCigarElements();
				for(int i = 0; i < (elList.size()-1); i++) {
					currentTel = elList.get(i);
					if(!(currentTel.getOperator().toString().equals("H") || currentTel.getOperator().toString().equals("D")))
						trackerSeq += elList.get(i).getLength();
				}
				for(int i = 0; i < tel.getLength(); i++) {
					if(seqRead[trackerSeq] == 'A')
						polyA++;
					else if (seqRead[trackerSeq] == 'T')
						polyT++;
					trackerSeq++;
				}
			}
			else
				continue;
			if((polyA / (double)tel.getLength()) >= 0.5 || (polyT/(double)tel.getLength()) >= 0.5)
				count++;
		}
		it.close();
		return ((count/total) >= 0.5);
	}
	private boolean rteEvent(svInfo current, svInfo next, File bam, File reference, int mini) {
		System.out.println("Entered RTE analysis...");
		SAMRecordIterator it = null;
		SamReader inputSam = null;
		if(bam.toString().endsWith(".cram") && CramIO.checkHeaderAndEOF(bam))
			inputSam = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).referenceSequence(reference).open(SamInputResource.of(bam));
		else
			inputSam = SamReaderFactory.make().enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX).validationStringency(ValidationStringency.LENIENT).samRecordFactory(DefaultSAMRecordFactory.getInstance()).open(bam);
		if(!intersect(current, next, mini)) {
			QueryInterval[] intervals = {new QueryInterval(inputSam.getFileHeader().getSequenceIndex(current.getChromLeft()), current.getStartLeft()-mini, current.getEndLeft() + mini), new QueryInterval(inputSam.getFileHeader().getSequenceIndex(next.getChromLeft()), next.getStartLeft()-mini, next.getEndLeft() + mini)};
			it = inputSam.query(intervals, false);
		} else
			it = inputSam.query(current.getChromLeft(), current.getStartLeft() - mini, next.getEndLeft() + mini, false);
		SAMRecord currentRecord = null;
		double count = 0;
		double total = 0;
		readloop:
		while(it.hasNext()) {
			double polyA = 0;
			double polyT = 0;
			boolean first = false;
			int trackerSeq = 0;
			currentRecord = it.next();
			byte[] seqRead = currentRecord.getReadBases();
			Iterator<CigarElement> cigIterator = currentRecord.getCigar().getCigarElements().iterator();
			CigarElement tel = null;
			try {
				while(cigIterator.hasNext()) {
					tel = cigIterator.next();
					if(tel.getOperator().toString().equals("H"))
						continue;
					else if(tel.getOperator().toString().equals("S")) {
						if(tel.getLength() < 20) {
							trackerSeq += tel.getLength();
							continue;
						} else {
							if(!first) {
								first = true;
								total++;
							}
							for(int i = 0; i < tel.getLength(); i++) {
								if(seqRead[trackerSeq] == 'A')
									polyA++;
								else if (seqRead[trackerSeq] == 'T')
									polyT++;
								trackerSeq++;
							}
						}
						if((polyA / (double)tel.getLength()) >= 0.5 || (polyT/(double)tel.getLength()) >= 0.5) {
							count++;
							continue readloop;
						}
					} else if(tel.getOperator().toString().equals("M")) 
						trackerSeq += tel.getLength();
					else if(tel.getOperator().toString().equals("I"))
						trackerSeq += tel.getLength();
					else if(tel.getOperator().toString().equals("D")) {}
					else {
						System.out.println("You are missing the following operator: " + tel.getOperator().toString());
						System.exit(0);
					}
				}
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(-2);
			}
		}
		return ((count/total) >= 0.5);
	}
	private boolean intersect(svInfo current, svInfo next, int search) {
		if(Math.abs(current.getStartLeft() - next.getStartLeft()) <= (2*search))
			return true;
		return false;
	}
	
	private ArrayList<svInfo> getStructuralVariants(File bed) {
		ArrayList<svInfo> svs = new ArrayList<svInfo>(10000);
		try {
			var br = new BufferedReader(new FileReader(bed));
			String line = null;
			String[] tokens = null;
			while((line = br.readLine()) != null) {
				if(line.startsWith("#") || line.contains("deletion") || line.contains("tandem"))
					continue;
				tokens = line.split("\t");
				svs.add(new svInfo(tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]), tokens[3], Integer.parseInt(tokens[4]), Integer.parseInt(tokens[5]), tokens[8].equals("+"), tokens[9].equals("+"), tokens[11], Integer.parseInt(tokens[12])));
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return svs;
	}
}