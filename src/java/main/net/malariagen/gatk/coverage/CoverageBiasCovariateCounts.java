package net.malariagen.gatk.coverage;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;

import net.malariagen.utils.io.TsvWriter;


public class CoverageBiasCovariateCounts implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 6408560452914950488L;
	public static final String SERIAL_FILE_NAME = "serialized.jso";
	public static final String RATES_FILE_NAME = "rates.tsv";
	public static final String SITE_COUNTS_FILE_NAME = "sites.tsv";
	public static final String STARTS_FILE_NAME = "starts.tsv";
	private Map<String,Map<String,CoverageBiasCovariateCounts.CovariateCombination>> locusCovariates;
	
	CoverageBiasCovariateCounts(Set<String> groupNames) {
		locusCovariates = new HashMap<String,Map<String,CoverageBiasCovariateCounts.CovariateCombination>>(groupNames.size());
		for (String g : groupNames)
			locusCovariates.put(g,new HashMap<String,CoverageBiasCovariateCounts.CovariateCombination>());
	}
	
	void add(String group, double gcBias, int size, long sites, long starts) {
		Map<String,CoverageBiasCovariateCounts.CovariateCombination> lcm = locusCovariates.get(group);
		if (lcm == null)
			throw new IllegalArgumentException("group " + group + " was not provided during contruction");
		String key = CovariateCombination.keyFor(gcBias, size);
		CoverageBiasCovariateCounts.CovariateCombination lc = lcm.get(key);
		if (lc == null)
			lcm.put(key, lc = new CovariateCombination(gcBias, size));
		lc.increment(sites,starts);
	}
	
	public void mergeIn(CoverageBiasCovariateCounts other) {
		if (other == null)
			throw new IllegalArgumentException("the other counts cannot be null");
		for (Map.Entry<String, Map<String, CovariateCombination>> e : other.locusCovariates.entrySet())
			if (!this.locusCovariates.containsKey(e.getKey()))
				throw new IllegalArgumentException("the other counts contain unknown group name " + e.getKey());
			else {
				Map<String,CovariateCombination> mineMap = locusCovariates.get(e.getKey());
				for (CovariateCombination cc : e.getValue().values()) {
					CovariateCombination mine = mineMap.get(cc.key());
					if (mine == null) {
						mineMap.put(cc.key(),cc.clone());
					}
					else {
						mine.siteCount += cc.siteCount;
						mine.startCount += cc.startCount;
					}
				}
			}
	}		
	
	public void saveIn(File outDir) throws IOException {
		outDir.mkdir();
		File serFile = new File(outDir,SERIAL_FILE_NAME);
		ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(serFile));
		oos.writeObject(this);
		oos.close();
		Set<String> groupNames = locusCovariates.keySet();
		TreeSet<CovariateCombination> allCC = new TreeSet<CovariateCombination>();
		for (Map<String,CovariateCombination> cc : this.locusCovariates.values())
			allCC.addAll(cc.values());
		
		String[] header = new String[groupNames.size() + 2];//+ 3];
		header[0] = "GCBias"; header[1] = "Size";
		File tableFile = new File(outDir,RATES_FILE_NAME);
		File siteCountsFile = new File(outDir,SITE_COUNTS_FILE_NAME);
		File startCountsFile = new File(outDir,STARTS_FILE_NAME);
		TsvWriter tw = new TsvWriter(new FileWriter(tableFile));
		TsvWriter scTw = new TsvWriter(new FileWriter(siteCountsFile));
		TsvWriter stTw = new TsvWriter(new FileWriter(startCountsFile));
		int nextIdx = 2; //= 3;
		for (String gn : groupNames)
			header[nextIdx++] = gn;
		tw.writeLine((Object[])header);
		scTw.writeLine((Object[])header);
		stTw.writeLine((Object[])header);
		Object[] values = new Object[header.length];
		Object[] siteCounts = new Object[header.length];
		Object[] startCounts = new Object[header.length];
		for (CovariateCombination ck : allCC) {
			String key = ck.key();
			values[0] = siteCounts[0] = startCounts[0] = String.format("%.2f",ck.gcBias);
 			values[1] = siteCounts[1] = startCounts[1] = "" + ck.size;
			for (int i = 2; i < header.length; i++) {
				CovariateCombination cc = locusCovariates.get(header[i]).get(key);
				values[i] = cc == null ? "NA" : String.format("%.4f",(double) cc.startCount / (double) cc.siteCount);
				siteCounts[i] = cc == null ? 0 : cc.siteCount;
				startCounts[i] = cc == null ? "NA" : cc.startCount; 
			}
			scTw.writeLine(siteCounts);
			stTw.writeLine(startCounts);
			tw.writeLine(values);
		}
		scTw.close();
		stTw.close();
		tw.close();
	}
	
	public static CoverageBiasCovariateCounts loadFrom(File outDir) throws IOException {
		if (outDir == null)
			throw new IllegalArgumentException("storage file or directory must not be null");
		if (!outDir.exists())
			throw new IllegalArgumentException("storage file or directory must exists (" + outDir + ")");
		File serializedFile = outDir.isDirectory() ? 
			serializedFile = new File(outDir,SERIAL_FILE_NAME) : outDir;
			
		ObjectInputStream ois = new ObjectInputStream(new FileInputStream(serializedFile));
		try {
			return (CoverageBiasCovariateCounts) ois.readObject();
		} catch (ClassNotFoundException e) {
			throw new IllegalArgumentException("the file provided does not contain coverage-bias covariate counts",e);
		}
	}

	
	static class CovariateCombination implements Serializable, Comparable<CovariateCombination>, Cloneable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		public final double gcBias;
		public long siteCount;
		public long startCount;
		private int hashCode = -1;
		public final int size;
		
		@Override
		public CovariateCombination clone() {
			try {
				return (CovariateCombination) super.clone();
			} catch (CloneNotSupportedException e) {
				throw new RuntimeException("this must never happen",e);
			}

		}
		
		@Override
		public boolean equals(Object o) {
			if (o == null)
				return false;
			if (!(o instanceof CovariateCombination))
				return false;
			CovariateCombination other = (CovariateCombination)o;
			return compareTo(other) == 0;
		}
		
		public String key() {
			return keyFor(gcBias, size);
		}

		@Override
		public int compareTo(CovariateCombination other) {
			int result = 0;
			if ((result = Double.compare(other.gcBias,this.gcBias)) != 0)
				return result;
			if ((result = Integer.compare(other.size,this.size)) != 0)
				return result;
			return 0;
		}
		
		@Override
		public int hashCode() {
			if (hashCode == -1)
				hashCode = keyFor(gcBias,size).hashCode();
			return hashCode;
		}
		
		public static String keyFor(double gcBias, int size) {
//			return String.format("%.2f %.2f %.2f",gcBias,nucEnt,trinucEnt);
			return String.format("%.2f %d",gcBias,size);
		}
		
		public void increment(long sites, long starts) {
			siteCount += sites;
			startCount += starts;
		}

		public CovariateCombination(double gcBias, int size) {
			this.gcBias = gcBias;
			this.size = size;
			siteCount = 0;
			startCount = 0;
		}
		

	}
	

}