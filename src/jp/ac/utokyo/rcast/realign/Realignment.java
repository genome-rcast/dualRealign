/*
 * Copyright Hiroki Ueda

 *  This program is free software; you can redistribute it and/or modify it under
 *	the terms of the GNU General Public License as published by the Free Software
 *	Foundation; either version 2 of the License, or (at your option) any later
 *	version.
	
 *	This program is distributed in the hope that it will be useful, but WITHOUT
 *	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *	FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 *	details.
	
 *	You should have received a copy of the GNU General Public License along with
 *	this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 *	Place, Suite 330, Boston, MA 02111-1307 USA
 */


package jp.ac.utokyo.rcast.realign;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.TumorGenotyper;
import jp.ac.utokyo.rcast.karkinos.utils.OptionComparator;
import jp.ac.utokyo.rcast.karkinos.utils.ReadWriteBase;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang3.StringUtils;

import srma.SRMA_mod;
import srma.TwobitReferenceSequence;

public class Realignment extends ReadWriteBase {

	public static void main(String[] arg) {

		//
		BasicParser parcer = new BasicParser();
		List<Option> optionList = getOptionListForKarkinos();
		Options opts = new Options();
		for (Option opt : optionList) {
			opts.addOption(opt);
		}

		CommandLine cl = null;
		try {
			cl = parcer.parse(opts, arg);
		} catch (ParseException e1) {
			System.out.println(e1.getMessage());
			HelpFormatter help = new HelpFormatter();
			help.setOptionComparator(new OptionComparator(optionList));
			help.printHelp("karkinos.jar analysis", opts, true);
			return;
		}

		Realignment ra = new Realignment();
		try {
			ra.exec(cl);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private void exec(CommandLine cl) throws Exception {

		String normalbamf = cl.getOptionValue("n");
		String tumorbamf = cl.getOptionValue("t");
		String twobitref = cl.getOptionValue("r");

		String targetRegion = null;
		if (cl.hasOption("ct")) {
			targetRegion = cl.getOptionValue("ct");
		}
		String outdir = cl.getOptionValue("o");
		if (!outdir.endsWith("/")) {
			outdir = outdir + "/";
			File f = new File(outdir);
			if (!f.exists()) {
				boolean suc = false;
				try {
					suc = f.mkdirs();
				} catch (Exception ex) {
					System.out.println("could not make directory " + outdir);
					return;
				}
				if (suc == false) {
					System.out.println("could not make directory " + outdir);
					return;
				}
			}
		}

		int numthread = 1;
		if (cl.hasOption("nt")) {
			numthread = Integer.parseInt(cl.getOptionValue("nt"));
		}

		try {
			realign(normalbamf, tumorbamf, twobitref, targetRegion, outdir,
					numthread);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private void realign(String normalbamf, String tumorbamf, String twobitref,
			String targetRegion, String outdir, int numthread) throws Exception {

		// load target
		TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(twobitref));
		DataSet dataset = null;
		if (targetRegion != null) {
			dataset = new DataSet(tgr.readIndex());
			dataset.loadTargetBed(targetRegion, tgr);
		}
		//
		
		SAMFileReader forheader = getReader(normalbamf);
		String outnameN = outdir + "/" + new File(normalbamf).getName()
				+ "_realign.bam";
		String outnameT = outdir + "/" + new File(tumorbamf).getName()
				+ "_realign.bam";

		SAMFileWriter normalbamw = getPreSortWriter(forheader.getFileHeader(),
				outnameN);
		SAMFileWriter tumorbamw = getPreSortWriter(forheader.getFileHeader(), outnameT);

		List<SAMSequenceRecord> ssrList = forheader.getFileHeader()
				.getSequenceDictionary().getSequences();
		forheader.close();
		
		
		int chromcnt = 0;

		System.out.println("start realgin");
		int cnt = 0;
		for (SAMSequenceRecord ssr : ssrList) {

			String chrom = ssr.getSequenceName();
			boolean usualchrom = usualChrom(chrom);

			System.out.println("processing chr " + chrom);
			
			SAMFileReader normalbamr = getReader(normalbamf);
			SAMFileReader normalbamr2 = getReader(normalbamf);
			SAMFileReader tumorbamr = getReader(tumorbamf);
			SAMFileReader tumorbamr2 = getReader(tumorbamf);
			if (usualchrom) {

				
				realignChrom(chrom, normalbamr, normalbamw, tumorbamr,
						tumorbamw, tgr, dataset, numthread,normalbamr2,tumorbamr2);

			} else {

				
				copy(chrom, normalbamr, normalbamw);
				copy(chrom, tumorbamr, tumorbamw);

			}
			normalbamr.close();
			tumorbamr.close();
			normalbamr2.close();
			tumorbamr2.close();

		}

		
		normalbamw.close();
		tumorbamw.close();

	}

	public static final int FlgNormal = 1;
	public static final int FlgTumor = 2;

	private void realignChrom(String chrom, SAMFileReader normalbamr,
			SAMFileWriter normalbamw, SAMFileReader tumorbamr,
			SAMFileWriter tumorbamw, TwoBitGenomeReader tgr, DataSet dataset,
			int numthread, SAMFileReader normalbamr2, SAMFileReader tumorbamr2) throws Exception {

				//
		TreeMap<Integer, Integer> indelpos = new TreeMap<Integer, Integer>();

		// stats indel pos
		System.out.println("stat indel pos " + chrom);

		CloseableIterator<SAMRecord> iteN = normalbamr
				.query(chrom, 0, 0, false);
		while (iteN.hasNext()) {

			SAMRecord sam = iteN.next();
			if (contatinIndel(sam)) {
				int ip = indelpos(sam);

				if (dataset != null) {
					if (null != dataset.getCh().getCapInterval(chrom, ip)) {

						if (indelpos.containsKey(ip)) {
							indelpos.put(ip, indelpos.get(ip) + 1);
						} else {
							indelpos.put(ip, 1);
						}

					}
				} else {

					if (indelpos.containsKey(ip)) {
						indelpos.put(ip, indelpos.get(ip) + 1);
					} else {
						indelpos.put(ip, 1);
					}

				}
			}

		}
		iteN.close();

		CloseableIterator<SAMRecord> iteT = tumorbamr.query(chrom, 0, 0, false);
		while (iteT.hasNext()) {

			SAMRecord sam = iteT.next();
			if (contatinIndel(sam)) {
				int ip = indelpos(sam);

				if (dataset != null) {
					if (null != dataset.getCh().getCapInterval(chrom, ip)) {

						if (indelpos.containsKey(ip)) {
							indelpos.put(ip, indelpos.get(ip) + 1);
						} else {
							indelpos.put(ip, 1);
						}

					}
				} else {

					if (indelpos.containsKey(ip)) {
						indelpos.put(ip, indelpos.get(ip) + 1);
					} else {
						indelpos.put(ip, 1);
					}

				}

			}

		}
		iteT.close();

		Set<Integer> single = new HashSet<Integer>();
		int cnt = 0;
		for (Entry<Integer, Integer> et : indelpos.entrySet()) {

			if (et.getValue() >= 400) {
				single.add(et.getKey());
			}
			if (et.getValue() <= 2) {
				single.add(et.getKey());
			}

		}
		for (int key : single) {
			indelpos.remove(key);
		}

		List<SAMRecord> normal = new ArrayList<SAMRecord>();
		List<SAMRecord> tumor = new ArrayList<SAMRecord>();
		//
		List<SAMRecord> realgin = new ArrayList<SAMRecord>();

		//
		// to Normal, Tumor, mixFor realgin
		System.out.println("extract realgin reads " + chrom);

		CloseableIterator<SAMRecord> iteN2 = normalbamr2.query(chrom, 0, 0,
				false);
		while (iteN2.hasNext()) {

			SAMRecord sam = iteN2.next();
			if (contatinIndel(sam) || nearIndel(sam, indelpos)) {

				sam.setAttribute("YY", FlgNormal);
				realgin.add(sam);

			} else {

				normal.add(sam);

			}

		}
		iteN2.close();

		CloseableIterator<SAMRecord> iteT2 = tumorbamr2
				.query(chrom, 0, 0, false);

		while (iteT2.hasNext()) {

			SAMRecord sam = iteT2.next();
			if (contatinIndel(sam) || nearIndel(sam, indelpos)) {

				sam.setAttribute("YY", FlgTumor);
				realgin.add(sam);

			} else {

				tumor.add(sam);

			}

		}
		iteT2.close();

		if (realgin.size() > 1) {

			System.out.println("normal not to realgin " + normal.size());
			System.out.println("tumor not to realgin " + tumor.size());

			System.out.println("sort reads for realgin " + chrom);
			// sort realgin
			Collections.sort(realgin, new SAMRecordCoordinateComparator());

			int start = realgin.get(0).getAlignmentStart();
			int end = realgin.get(realgin.size() - 1).getAlignmentEnd();
			//

			System.out.println("realign by SRMA " + chrom + " size ="
					+ realgin.size());

			TwobitReferenceSequence tbrs = new TwobitReferenceSequence(tgr,
					normalbamr.getFileHeader());
			ReferenceSequence res = tbrs.getSequence(chrom, start, end);

			// realgin
			// ////////////////////////////////////

			// sep list
			List<List<SAMRecord>> sep = sep(realgin);
			List<List<SAMRecord>> ret = new ArrayList<List<SAMRecord>>();
			ExecutorService exec = Executors.newFixedThreadPool(numthread);

			try {

				for (List<SAMRecord> list : sep) {

					Runtask task = new Runtask(list, ret,
							normalbamr.getFileHeader(), res, tbrs, chrom,
							start, end);

					exec.execute(task);

				}

			} finally {
				exec.shutdown();
				exec.awaitTermination(5, TimeUnit.HOURS);
			}

			System.out.println("realign by SRMA finish " + chrom);

			// ///////////////////////////////////
			//
			for (List<SAMRecord> rlist : ret) {
				// add to each
				for (SAMRecord sam : rlist) {

					int flg = sam.getIntegerAttribute("YY");
					if (flg == FlgNormal) {
						normal.add(sam);
					} else {
						tumor.add(sam);
					}

				}
			}

			System.out.println("sort normal " + chrom);
			// sort normal
			Collections.sort(normal, new SAMRecordCoordinateComparator());
			for (SAMRecord sam : normal) {

				normalbamw.addAlignment(sam);

			}
			normal = null;
			// normalbamw.close();

			System.out.println("sort tumor " + chrom);
			Collections.sort(tumor, new SAMRecordCoordinateComparator());
			for (SAMRecord sam : tumor) {

				tumorbamw.addAlignment(sam);

			}
			tumor = null;
			// tumorbamw.close();

		}

	}

	class Runtask implements Runnable {

		List<SAMRecord> list;
		List<List<SAMRecord>> ret;
		SAMFileHeader fileHeader;
		ReferenceSequence rsf;
		ReferenceSequenceFile referenceSequenceFile;
		String chr;
		int start;
		int end;

		Runtask(List<SAMRecord> list, List<List<SAMRecord>> ret,
				SAMFileHeader fileHeader, ReferenceSequence rsf,
				ReferenceSequenceFile referenceSequenceFile, String chr,
				int start, int end) {

			this.list = list;
			this.ret = ret;
			this.fileHeader = fileHeader;
			this.rsf = rsf;
			this.referenceSequenceFile = referenceSequenceFile;
			this.chr = chr;
			this.start = start;
			this.end = end;

		}

		@Override
		public void run() {

			SRMA_mod inst = new SRMA_mod();
			try {
				List<SAMRecord> rlist = inst.doWorkSRMA(fileHeader, rsf,
						referenceSequenceFile, chr, start, end, list, 1);

				synchronized (ret) {
					ret.add(rlist);
				}

			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}

	public static int gap = 10000;

	private List<List<SAMRecord>> sep(List<SAMRecord> realgin) {

		List<List<SAMRecord>> l = new ArrayList<List<SAMRecord>>();

		int b4 = 0;
		List<SAMRecord> holder = new ArrayList<SAMRecord>();
		for (SAMRecord sam : realgin) {

			if (b4 == 0) {
				b4 = sam.getAlignmentEnd();
				continue;
			}

			if ((sam.getAlignmentStart() - b4) > gap) {

				//
				l.add(holder);
				holder = new ArrayList<SAMRecord>();
				holder.add(sam);

			} else {

				holder.add(sam);

			}

			b4 = sam.getAlignmentEnd();

		}
		l.add(holder);

		return l;
	}

	private boolean nearIndel(SAMRecord sam, TreeMap<Integer, Integer> indelpos) {

		int s = sam.getAlignmentStart();
		int e = sam.getAlignmentEnd();
		//
		int readlen = sam.getReadLength();

		Integer f = indelpos.ceilingKey(s);
		Integer c = indelpos.floorKey(e);

		if (f != null) {
			if (Math.abs(s - f) < readlen) {
				return true;
			}
		}
		if (c != null) {
			if (Math.abs(c - e) < readlen) {
				return true;
			}
		}
		return false;
	}

	private Integer indelpos(SAMRecord sam) {

		int pos = sam.getAlignmentStart();
		for (CigarElement ce : sam.getCigar().getCigarElements()) {

			if (ce.getOperator().equals(CigarOperator.D)) {
				return pos;
			}
			if (ce.getOperator().equals(CigarOperator.I)) {
				return pos;
			}

			if (ce.getOperator().consumesReferenceBases()) {
				pos = pos + ce.getLength();
			}
		}

		return pos;

	}

	private boolean contatinIndel(SAMRecord sam) {

		//
		//

		for (CigarElement ce : sam.getCigar().getCigarElements()) {

			if (ce.getOperator().equals(CigarOperator.D)) {
				return true;
			}
			if (ce.getOperator().equals(CigarOperator.I)) {
				return true;
			}

		}

		return false;

	}

	private void copy(String chrom, SAMFileReader bamr, SAMFileWriter bamw) {

		CloseableIterator<SAMRecord> ite = bamr.query(chrom, 0, 0, false);
		while (ite.hasNext()) {

			SAMRecord sam = ite.next();
			bamw.addAlignment(sam);

		}

	}

	private boolean usualChrom(String chrom) {

		String chromnum = chrom;
		if (chromnum.contains("chr")) {
			chromnum = chromnum.replaceAll("chr", "");
		}
		if (StringUtils.isNumeric(chromnum))
			return true;
		if (chromnum.equalsIgnoreCase("X"))
			return true;
		if (chromnum.equalsIgnoreCase("Y"))
			return true;

		return false;

	}

	public static Option getOption(String opt, String longOpt, boolean hasArg,
			String description, boolean required) {
		Option option = new Option(opt, longOpt, hasArg, description);
		option.setRequired(required);
		return option;
	}

	private static List<Option> getOptionListForKarkinos() {

		List<Option> optionlist = new ArrayList<Option>();
		optionlist.add(getOption("n", "normalBam", true, "normal bam file",
				true));
		optionlist
				.add(getOption("t", "tumorBam", true, "tumor bam file", true));

		optionlist.add(getOption("r", "reference", true,
				"2 bit genome reference file", true));

		optionlist
				.add(getOption("o", "outdir", true, "output directory", true));

		optionlist.add(getOption("ct", "captureTarget", true,
				"Capture target regions(bed format)", false));

		optionlist.add(getOption("nt", "num threads", true, "number of threads",
				false));

		return optionlist;

	}

}
