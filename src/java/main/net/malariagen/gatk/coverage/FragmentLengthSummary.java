package net.malariagen.gatk.coverage;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import net.malariagen.gatk.math.EmptyIntegerDistribution;
import net.malariagen.gatk.math.IntegerDistribution;


public class FragmentLengthSummary {

	Map<String,IntegerDistribution> smFlDist;
	Map<String,IntegerDistribution> smIlDist;
	Map<String,IntegerDistribution> rgFlDist;
	Map<String,IntegerDistribution> rgIlDist;
	
	
	FragmentLengthSummary(FragmentLengthFrequencies freqs) {
		this.smFlDist = new LinkedHashMap<String,IntegerDistribution>(freqs.smIndex.size());
		this.smIlDist = new LinkedHashMap<String,IntegerDistribution>(freqs.smIndex.size());
		this.rgFlDist = new LinkedHashMap<String,IntegerDistribution>(freqs.rgIndex.size());
		this.rgIlDist = new LinkedHashMap<String,IntegerDistribution>(freqs.rgIndex.size());
		
		buildSampleDistributions(freqs);
		buildRgDistributions(freqs);
			
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
				if (rgFlFreqs[i] != 0) diffCount++;
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

			IntegerDistribution rgFlDist = values.length == 0 ? new EmptyIntegerDistribution() : new MyIntegerDistribution(values,accu);
			this.rgFlDist.put(s, rgFlDist);
			long[] rgIlFreqs = freqs.rgIlFreq[rgIdx];
			diffCount = 0;
			for (int i = 1; i < rgIlFreqs.length; i++)
				if (rgIlFreqs[i] != 0) diffCount++;
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
			IntegerDistribution rgIlDist = values.length == 0 ? new EmptyIntegerDistribution() : new MyIntegerDistribution(values,accu);
			this.rgIlDist.put(s,rgIlDist);
		}	
	}

	private void buildSampleDistributions(FragmentLengthFrequencies freqs) {
		for (String s : freqs.samples) {
			int smIdx = freqs.smIndex.get(s);
			long[] smFlFreqs = freqs.smFlFreq[smIdx];
			int diffCount = 0;
			for (int i = 1; i < smFlFreqs.length; i++)
				if (smFlFreqs[i] != 0) diffCount++;
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
			
			IntegerDistribution smFlDist = values.length == 0 ? new EmptyIntegerDistribution() : new MyIntegerDistribution(values,accu);
			this.smFlDist.put(s, smFlDist);
			
			long[] smIlFreqs = freqs.smIlFreq[smIdx];
			diffCount = 0;
			for (int i = 1; i < smIlFreqs.length; i++)
				if (smIlFreqs[i] != 0) diffCount++;
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
			IntegerDistribution smIlDist = values.length == 0 ? new EmptyIntegerDistribution() : new MyIntegerDistribution(values,accu);
			this.smIlDist.put(s,smIlDist);
		}
	}
	
	public class MyIntegerDistribution implements IntegerDistribution {

		int[] values;
		long[] accu;
		int mean = -1;
		int var = -1;
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
			if (mean == -1) calculateMean();
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
			long offset = Math.round(q * count);
			int where = Arrays.binarySearch(accu, offset);
			if (where < 0) {
				where = - (where + 1);
			}
			else if (where == accu.length) {
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
			return values[0];
		}

		@Override
		public int maximum() {
			return values[values.length - 1];
		}

		@Override
		public double cumulativeProbability(int c) {
			int where = Arrays.binarySearch(values, c);
			if (where < 0) {
				where = - (where + 1);
			}
			else if (where == accu.length) {
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
			int bestMode = 0;
			long sum = 0;
			double ssum = 0;
			for (int i = 0; i < accu.length; i++) {
				long a = accu[i];
				int v = values[i];
				if (a > accu[bestMode])
					bestMode = i;
				sum += a * v;
				ssum += Math.log(a * v * a * v);
			}
			mean = (int) Math.round((double) sum / (double) count());
			var = (int) Math.abs(mean * mean - ssum / (double) count());
			mode = values[bestMode];
		}
		
	}
}
