package net.malariagen.gatk.coverage;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Reader;
import java.io.Serializable;


import java.util.Arrays;

import java.util.LinkedHashMap;
import java.util.Map;


import net.malariagen.gatk.math.EmptyIntegerDistribution;
import net.malariagen.gatk.math.IntegerDistribution;
import net.malariagen.utils.io.TsvWriter;

public class FragmentLengthSummary implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 4519497188985297029L;

	public static final String SUMMARY_FILE_NAME = "summary";
	
	Map<String, MyIntegerDistribution> smFlDist;
	Map<String, MyIntegerDistribution> smIlDist;
	Map<String, MyIntegerDistribution> rgFlDist;
	Map<String, MyIntegerDistribution> rgIlDist;

	FragmentLengthSummary(FragmentLengthFrequencies freqs) {
		this.smFlDist = new LinkedHashMap<String, MyIntegerDistribution>(
				freqs.smIndex.size());
		this.smIlDist = new LinkedHashMap<String, MyIntegerDistribution>(
				freqs.smIndex.size());
		this.rgFlDist = new LinkedHashMap<String, MyIntegerDistribution>(
				freqs.rgIndex.size());
		this.rgIlDist = new LinkedHashMap<String, MyIntegerDistribution>(
				freqs.rgIndex.size());

		buildSampleDistributions(freqs);
		buildRgDistributions(freqs);

	}

	public static FragmentLengthSummary loadFrom(File  dirIn) throws IOException {
		if (dirIn == null) 
			throw new IllegalArgumentException("the input directory cannot be null");
		if (!dirIn.exists())
			throw new IOException("the input directory does not exists or is not reachable");
		if (!dirIn.isDirectory())
			throw new IOException("the input directory is not in fact a a directory");
		if (!dirIn.canRead())
			throw new IOException("the input directory does not allow to read from");
		
		File serializeFile = new File(dirIn,"serialized.jso");
		if (!serializeFile.exists())
			throw new IOException("the input serialize file does not exists or is not reachable");
		if (!serializeFile.isFile())
			throw new IOException("the input serialize file is not a regular file");
		if (!serializeFile.canRead())
			throw new IOException("the input serialize file cannot be read");
		
		ObjectInputStream ois = new ObjectInputStream(new FileInputStream(serializeFile));
		Object result;
		try {
			result = ois.readObject();
		} catch (ClassNotFoundException e) {
			throw new IOException("serialize file content is not a fragment length summary or of the wrong version",e);
		}
		if (!(result instanceof FragmentLengthSummary)) 
			throw new IOException("serialize file content is not a fragment length summary");
		return (FragmentLengthSummary) result;
	}

	public void saveIn(File dirOut) throws IOException {
		File summaryFile = new File(dirOut, SUMMARY_FILE_NAME);
		File histogramFile = new File(dirOut, "histogram.tsv");
		File serializedFile = new File(dirOut, "serialized.jso");
		dirOut.mkdir();
		if (!dirOut.exists())
			throw new IOException(
					"could not create output directory for fragment length summary files");
		else if (!dirOut.isDirectory())
			throw new IOException("existing file '" + dirOut
					+ " is not a directory");
		else if (!dirOut.canWrite())
			throw new IOException(
					"cannot write in output directory for fragment length summary files");
		saveSummaryFileIn(summaryFile);
		saveHistogramFileIn(histogramFile);
		saveSerializedIn(serializedFile);
	}

	private void saveSerializedIn(File serializedFile) throws IOException {
		ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(serializedFile));
		oos.writeObject(this);
		oos.close();
	}

	private void saveHistogramFileIn(File histogramFile) throws IOException {
		TsvWriter tsv = new TsvWriter(new FileWriter(histogramFile));
		int rgCount = rgFlDist.size();
		int smCount = smFlDist.size();
		int allCount = rgCount + smCount + 1;
		Object[] header = new String[allCount];
		header[0] = "LENGTH";
		MyIntegerDistribution[] integerDistributions = new MyIntegerDistribution[allCount - 1];
		int nextIdx = 1;
		int maxLength = Integer.MIN_VALUE;
		int minLength = Integer.MAX_VALUE;
		for (Map.Entry<String,MyIntegerDistribution> e : rgFlDist.entrySet()) {
			MyIntegerDistribution dist = e.getValue();
			if ((integerDistributions[nextIdx-1] = dist).maximum() > maxLength)
				maxLength = dist.maximum();
			if (dist.minimum() < minLength)
				minLength = dist.minimum();
			header[nextIdx++] = e.getKey();
		}
		for (Map.Entry<String,MyIntegerDistribution> e : smFlDist.entrySet()) {
			MyIntegerDistribution dist = e.getValue();
			if ((integerDistributions[nextIdx-1] = dist).maximum() > maxLength)
				maxLength = dist.maximum();
			if (dist.minimum() < minLength)
				minLength = dist.minimum();
			header[nextIdx++] = e.getKey();
		}
		tsv.writeLine(header);
		int[] indeces = new int[allCount - 1];
		long[] values = new long[allCount];
		for (int len = minLength; len <= maxLength; len++) {
			values[0] = len;
			for (int d = 0; d < allCount - 1; d++) {
				int[] v = integerDistributions[d].values;
				long[] a = integerDistributions[d].accu;
				while (v.length > indeces[d] && v[indeces[d]] < len) indeces[d]++;
				if (v.length < indeces[d])
					continue;
				if (v[indeces[d]] == len)
					values[1+d] = a[indeces[d]] - (indeces[d] == 0 ? 0 : a[indeces[d] - 1]);
				else
					values[1+d] = 0;
			}
			tsv.writeLine(values);
		}
		tsv.close();
	}

	private void saveSummaryFileIn(File summaryFile) throws IOException {
		TsvWriter w = new TsvWriter(new FileWriter(summaryFile));
		w.writeLine("ID","TYPE","COUNT","MEAN","VAR","MIN","Q25","Q50","Q75","MAX");
		for (Map.Entry<String,MyIntegerDistribution> e : rgFlDist.entrySet()) {
			MyIntegerDistribution dist = e.getValue();
			w.writeLine(e.getKey(),"RG",dist.count(),dist.mean(),dist.variance(),dist.minimum(),dist.quantile(0.25),
					dist.quantile(0.50),dist.quantile(0.75),dist.maximum());
		}
		for (Map.Entry<String,MyIntegerDistribution> e : smFlDist.entrySet()) {
			MyIntegerDistribution dist = e.getValue();
			w.writeLine(e.getKey(),"SM",dist.count(),dist.mean(),dist.variance(),dist.minimum(),dist.quantile(0.25),
					dist.quantile(0.50),dist.quantile(0.75),dist.maximum());
		}
		w.close();
		
		
	}

	public IntegerDistribution getSampleInsertLengths(String s) {
		return smIlDist.get(s);
	}

	public IntegerDistribution getSampleFragmentLengths(String s) {
		return smFlDist.get(s);
	}

	public IntegerDistribution getReadGroupInsertLengths(String s) {
		return rgIlDist.get(s);
	}

	public IntegerDistribution getReadGroupFragmentLengths(String s) {
		return rgFlDist.get(s);
	}

	private void buildRgDistributions(FragmentLengthFrequencies freqs) {
		for (String s : freqs.readGroups) {
			int rgIdx = freqs.rgIndex.get(s);
			long[] rgFlFreqs = freqs.rgFlFreq[rgIdx];
			int diffCount = 0;
			for (int i = 1; i < rgFlFreqs.length; i++)
				if (rgFlFreqs[i] != 0)
					diffCount++;
			int[] values = new int[diffCount];
			long[] accu = new long[diffCount];
			int nextIdx = 0;
			long soFar = 0;
			for (int i = 1; i < rgFlFreqs.length; i++) {
				if (rgFlFreqs[i] == 0)
					continue;
				values[nextIdx] = i;
				soFar = accu[nextIdx++] = soFar + rgFlFreqs[i];
			}

			MyIntegerDistribution rgFlDist = new MyIntegerDistribution(values, accu);
			this.rgFlDist.put(s, rgFlDist);
			long[] rgIlFreqs = freqs.rgIlFreq[rgIdx];
			diffCount = 0;
			for (int i = 1; i < rgIlFreqs.length; i++)
				if (rgIlFreqs[i] != 0)
					diffCount++;
			values = new int[diffCount];
			accu = new long[diffCount];
			nextIdx = 0;
			soFar = 0;
			for (int i = 1; i < rgIlFreqs.length; i++) {
				if (rgIlFreqs[i] == 0)
					continue;
				values[nextIdx] = i;
				soFar = accu[nextIdx++] = soFar + rgIlFreqs[i];
			}
			MyIntegerDistribution rgIlDist = new MyIntegerDistribution(values, accu);
			this.rgIlDist.put(s, rgIlDist);
		}
	}

	private void buildSampleDistributions(FragmentLengthFrequencies freqs) {
		for (String s : freqs.samples) {
			int smIdx = freqs.smIndex.get(s);
			long[] smFlFreqs = freqs.smFlFreq[smIdx];
			int diffCount = 0;
			for (int i = 1; i < smFlFreqs.length; i++)
				if (smFlFreqs[i] != 0)
					diffCount++;
			int[] values = new int[diffCount];
			long[] accu = new long[diffCount];
			int nextIdx = 0;
			long soFar = 0;
			for (int i = 1; i < smFlFreqs.length; i++) {
				if (smFlFreqs[i] == 0)
					continue;
				values[nextIdx] = i;
				soFar = accu[nextIdx++] = soFar + smFlFreqs[i];
			}

			MyIntegerDistribution smFlDist = new MyIntegerDistribution(values, accu);
			this.smFlDist.put(s, smFlDist);

			long[] smIlFreqs = freqs.smIlFreq[smIdx];
			diffCount = 0;
			for (int i = 1; i < smIlFreqs.length; i++)
				if (smIlFreqs[i] != 0)
					diffCount++;
			values = new int[diffCount];
			accu = new long[diffCount];
			nextIdx = 0;
			soFar = 0;
			for (int i = 1; i < smIlFreqs.length; i++) {
				if (smIlFreqs[i] == 0)
					continue;
				values[nextIdx] = i;
				soFar = accu[nextIdx++] = soFar + smIlFreqs[i];
			}
			MyIntegerDistribution smIlDist = new MyIntegerDistribution(values, accu);
			this.smIlDist.put(s, smIlDist);
		}
	}

	public class MyIntegerDistribution implements IntegerDistribution, Serializable {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		int[] values;
		long[] accu;
		double mean = -1;
		double var = -1;
		int mode = -1;

		MyIntegerDistribution(int[] values, long[] accu) {
			this.values = values;
			this.accu = accu;
		}

		@Override
		public long count() {
			return this.accu[this.accu.length - 1];
		}

		@Override
		public double mean() {
			if (mean == -1)
				calculateMean();
			return mean;
		}

		private void calculateMean() {
			calculateMeanAndVar();
		}

		private void calculateVar() {
			calculateMeanAndVar();
		}

		private void calculateMode() {
			calculateMeanAndVar();
		}

		@Override
		public double variance() {
			if (var == -1)
				calculateVar();
			return var;
		}

		@Override
		public double standardDeviation() {
			return Math.sqrt(variance());
		}

		@Override
		public double percentile(int pc) {
			return quantile(pc * 0.01);
		}

		@Override
		public double quantile(double q) {
			long count = count();
			if (count == 0)
				return Double.NaN;
			long offset = Math.round(q * count);
			int where = Arrays.binarySearch(accu, offset);
			if (where < 0) {
				where = -(where + 1);
			} else if (where == accu.length) {
				where = accu.length - 1;
			}
			return values[where];
		}

		@Override
		public double median() {
			return quantile(0.5);
		}

		@Override
		public int minimum() {
			if (values.length == 0)
				return Integer.MIN_VALUE;
			return values[0];
		}

		@Override
		public int maximum() {
			if (values.length == 0)
				return Integer.MAX_VALUE;
			return values[values.length - 1];
		}

		@Override
		public double cumulativeProbability(int c) {
			if (values.length == 0)
				return Double.NaN;
			int where = Arrays.binarySearch(values, c);
			if (where < 0) {
				where = -(where + 1);
			} else if (where == accu.length) {
				where = accu.length - 1;
			}
			return (double) accu[where] / (double) count();
		}

		@Override
		public int mode() {
			if (mode == -1)
				calculateMode();
			return mode;
		}

		public void calculateMeanAndVar() {
			if (values.length == 0) {
				mean = Double.NaN;
				mode = 0;
				return;
			}
			int bestMode = 0;
			long sum = 0;
			double ssum = 0;
			long count = 0;
			for (int i = 0; i < accu.length; i++) {
				long a = accu[i];
				a -= count;
				int v = values[i];
				if (a > accu[bestMode])
					bestMode = i;
				count += a;
				sum += a * v;
				ssum += Math.log(a * v * a * v);
			}
			mean = (int) Math.round((double) sum / (double) count);
			var = (int) Math.abs(mean * mean - ssum / (double) count);
			mode = values[bestMode];
		}

	}

}
