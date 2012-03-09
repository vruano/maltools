package net.malariagen.gatk.math;

import java.util.Arrays;

import org.broadinstitute.sting.utils.collections.Pair;

import net.malariagen.gatk.coverage.AbstractCoverageDistribution;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonPrimitive;

class ValueFrequencyArrayIntegerDistribution extends
		AbstractCoverageDistribution {

	private int[] values;
	private long[] frequencies;
	private long count;
	private int modeIndex;
	private double mean;
	private double variance;
	private long sum;
	private long sqSum;

	private double[] cumulativeProbabilities;
	private double[] percentiles;
	// index in values of the value that correspond to that percentile.
	private int[] percentileIndices;
	// in what occurrence of the percentile value that percentile actually
	// falls.
	private long[] percentileFrequencyOffsets;
	
	ValueFrequencyArrayIntegerDistribution(long count, long sum, long sqSum,
			int[] values, long[] frequencies) {
		this.sum = sum;
		this.sqSum = sqSum;
		this.values = values;
		this.frequencies = frequencies;
		if (values.length != frequencies.length)
		  throw new IllegalArgumentException("value and frequency arrays must have the same lengths");
		this.count = count;
		mean = sum / (double) count;
		variance = sqSum / (double) count - (mean * mean);
		if (count == 0) {
			zeroCountSetup();
		} else {
			calculatePercentiles();
			mean = sum / (double) count;
			variance = sqSum / (double) count - mean * mean;
		}
	}

	private void zeroCountSetup() {
		count = 0;
		mean = variance = Double.NaN;
		modeIndex =  -1;
		values = new int[0];
		frequencies = new long[0];
		percentiles = new double[0];
		percentileIndices = new int[101];
		percentileFrequencyOffsets = new long[101];
		cumulativeProbabilities = new double[0];
	}

	// provided to be able to recover distributions from json files
	private ValueFrequencyArrayIntegerDistribution() {
	}

	private void calculatePercentiles() {

		PercentileSearchPointer psp = new PercentileSearchPointer();
		percentiles = new double[101]; // 0 to 100 inclusive so 101.
		percentileIndices = new int[101];
		percentileFrequencyOffsets = new long[101];
		cumulativeProbabilities = new double[values.length];
		for (int pc = 0; pc < 101; pc++) {
			percentiles[pc] = psp.seek(pc);
			percentileIndices[pc] = psp.valueIndex;
			percentileFrequencyOffsets[pc] = psp.frequencyOffset;
		}
	}

	class PercentileSearchPointer {
		int valueIndex = 0;
		long frequencyOffset = 0;
		long leftCount = 0;

		double seek(int pc) {
			if (pc < 0 || pc > 100)
				throw new IllegalArgumentException("percentage out of rage: " + pc);
			double pcFraction = pc == 100 ? 1.0 : ((double)pc) / 100.0;
			double exactBreak = count * pcFraction;
			long floorBreak = (long) Math.floor(exactBreak);
			long ceilBreak = floorBreak == count ? floorBreak : (long) Math.ceil(exactBreak);
			if (floorBreak < leftCount) 
			    reset();
			int floorValue = advance(floorBreak - leftCount);
			if (floorBreak == ceilBreak)
				return floorValue;
			int ceilValue = advance(ceilBreak - leftCount);
			return floorValue * (1 - (exactBreak - floorBreak)) + ceilValue
					* (1 - (ceilBreak - exactBreak));
		}

		void reset() {
			valueIndex = 0;
			frequencyOffset = 0;
			leftCount = 0;
		}

		// has a lot of side effects for efficiency (calculate cumulative-props. and mode) so please
		// modify with care.
		int advance(long increment) {
			while (increment > 0) {
				if (valueIndex >= frequencies.length)
					throw new IllegalStateException(" " + valueIndex + " " + frequencies.length + " " + leftCount + " " + count + " " + increment);
				long f = frequencies[valueIndex];
				long frequencyOffsetLeft = f - frequencyOffset;
				if (frequencyOffsetLeft > increment) {
					frequencyOffset += increment;
					leftCount += increment;
					increment = 0;
				}
				else {
					if (frequencies[modeIndex] < frequencies[valueIndex]) modeIndex = valueIndex;
					increment -= frequencyOffsetLeft;
					frequencyOffset = 0;
					leftCount += frequencyOffsetLeft;
					cumulativeProbabilities[valueIndex++] = ((double) leftCount)
							/ (double) count;
				}
			}
			
			int result = valueIndex >= values.length ? values[values.length - 1] : values[valueIndex];
			return result;
		    
		 }
	}

	@Override
	public double mean() {
		return mean;
	}

	@Override
	public double variance() {
		return variance;
	}

	@Override
	public double percentile(int pc) {
		return percentiles[pc];
	}

	@Override
	public double quantile(double q) {
		if (q < 0 || q > 1)
			throw new IllegalArgumentException("q cannot must be within [0,1]");
		int floor = (int) Math.floor(q * 100);
		int ceil = (int) Math.ceil(q * 100);
		if (floor == ceil) return percentiles[floor];
		long floorBreak = (long) Math.floor(q * count);
		long ceilBreak = (long) Math.ceil(q * count);
		PercentileSearchPointer psp = new PercentileSearchPointer();
		psp.valueIndex = percentileIndices[floor];
		psp.frequencyOffset = percentileFrequencyOffsets[floor];
		psp.leftCount = (long) Math.round(count * floor / 100);
		int floorValue = psp.advance(floorBreak -  psp.leftCount);
		if (floorBreak == ceilBreak)
			return floorValue;
		int ceilValue = psp.advance(ceilBreak - psp.leftCount);
		return floorValue * (1 - (q * count - floorBreak)) + ceilValue * (1 - (ceilBreak - q * count));
	}

	@Override
	public int minimum() {
		if (count == 0)
			throw new UnsupportedOperationException(
					"there is no minimum in a sample of zero elements");
		return values[0];
	}

	@Override
	public int maximum() {
		if (count == 0)
			throw new UnsupportedOperationException(
					"there is no maximum in a sample of zero elements");
		return values[values.length - 1];
	}

	@Override
	public double cumulativeProbability(int c) {
		if (c < minimum())
			return 0;
		else if (c > maximum())
			return 1;
		else {
			int index = Arrays.binarySearch(values, c);
			if (index >= 0)
				return cumulativeProbabilities[index];
			else {
				return cumulativeProbabilities[-index - 1];
			}
		}
	}
	
	public static ValueFrequencyArrayIntegerDistribution merge(ValueFrequencyArrayIntegerDistribution ... distributions) {
		if (distributions.length == 0)
			throw new IllegalArgumentException();
		if (distributions.length == 1)
			return distributions[0];
			
		Pair<int[],long[]> r = mergeArrays(distributions,0,distributions.length);
		long count = 0;
		long sum = 0;
		long sqSum = 0;
		for (ValueFrequencyArrayIntegerDistribution d : distributions) {
			count += d.count;
			sum += d.sum;
			sqSum += d.sqSum;
		}
		
		return new ValueFrequencyArrayIntegerDistribution(count,sum,sqSum,r.getFirst(),r.getSecond());
	}
	
	static Pair<int[],long[]> mergeArrays(ValueFrequencyArrayIntegerDistribution[] distributions, int offset, int length) {
		if (length == 0)
			return new Pair<int[],long[]>( new int[0], new long[0]);
		if (length == 1)
			return new Pair<int[],long[]>( distributions[offset].values , distributions[offset].frequencies );
		int length1 = length >> 1;
		int length2 = length - length1;
		Pair<int[],long[]> d1 = mergeArrays(distributions,offset,length1);
		Pair<int[],long[]> d2 = mergeArrays(distributions,offset + length1,length2);
		return mergeArrays(d1,d2);
	}
	
	static Pair<int[],long[]> mergeArrays(Pair<int[],long[]> d1, Pair<int[],long[]> d2) {
		int[] d1Values = d1.getFirst();
		long[] d1Frequencies = d1.getSecond();
		int[] d2Values = d2.getFirst();
		long[] d2Frequencies = d2.getSecond();
		int d1Length = d1Values.length;
		int d2Length = d2Values.length;		
		int[] values = new int[d1Length + d2Length];
		long[] frequencies = new long[d1Length + d2Length];
		
		int nextIndex = 0;
		int nextIndex1 = 0;
		int nextIndex2 = 0;
		
		while (nextIndex1 < d1Length && nextIndex2 < d2Length) {
			int v1 = d1Values[nextIndex1];
			int v2 = d2Values[nextIndex2];
			if (v1 < v2) {
				values[nextIndex] = v1;
				frequencies[nextIndex++] = d1Frequencies[nextIndex1++];
			}
			else if (v2 < v1) {
				values[nextIndex] = v2;
				frequencies[nextIndex++] = d2Frequencies[nextIndex2++];				
			}
			else { // same value.
				values[nextIndex] = v1; // == v2
				frequencies[nextIndex++] = d1Frequencies[nextIndex1++] + d2Frequencies[nextIndex2++];
			}
		}
		
		// finish off copying the trailing values for the still remaining input distribution if any.
		if (nextIndex1 < d1Length) {
			System.arraycopy(d1Values, nextIndex1, values, nextIndex, d1Length - nextIndex1);
			System.arraycopy(d1Frequencies, nextIndex1, frequencies, nextIndex, d1Length - nextIndex1);
		    nextIndex += d1Length - nextIndex1;
		}
		else if (nextIndex2 < d2Length) { 
			System.arraycopy(d2Values, nextIndex2, values, nextIndex, d2Length - nextIndex2);
			System.arraycopy(d2Frequencies, nextIndex2, frequencies, nextIndex, d2Length - nextIndex2);
		    nextIndex += d2Length - nextIndex2;
		}
		
		// truncate exceeding positions
		if (values.length != nextIndex) {
			values = Arrays.copyOf(values, nextIndex);
			frequencies = Arrays.copyOf(frequencies, nextIndex);
		}
		return new Pair<int[],long[]>(values, frequencies);
	}
	
	JsonObject toJsonObject() {
		JsonObject result = new JsonObject();
		result.add("count", new JsonPrimitive(count));
		if (count == 0)
			return result;
		result.add("mean", new JsonPrimitive(mean));
		result.add("mode", new JsonPrimitive(values[modeIndex]));
		result.add("variance", new JsonPrimitive(variance));
		result.add("sum", new JsonPrimitive(sum));
		result.add("sqSum", new JsonPrimitive(sqSum));
		result.add("values",toJsonArray(values));
		result.add("frequencies", toJsonArray(frequencies));
		result.add("cumulativeProbabilities", toJsonArray(cumulativeProbabilities));
		result.add("percentiles", toJsonArray(percentiles));
		result.add("percentileIndices", toJsonArray(percentileIndices));
		result.add("percentileFrequencyOffsets", toJsonArray(percentileFrequencyOffsets));
		return result;
		
	}
	
	
	static ValueFrequencyArrayIntegerDistribution fromJsonObject(JsonObject json) {
		ValueFrequencyArrayIntegerDistribution result = new ValueFrequencyArrayIntegerDistribution();
		result.loadJsonObject(json);
		return result;
	}
	
	private void loadJsonObject(JsonObject json) {
		count = json.get("count").getAsInt();
		if (count == 0) {
		  zeroCountSetup();
		}
		else {
		  mean = json.get("mean").getAsDouble();
		  variance = json.get("variance").getAsDouble();
		  sum = json.get("sum").getAsLong();
		  sqSum = json.get("sqSum").getAsLong();
		  values = fromJsonIntegerArray(json.get("values").getAsJsonArray());
		  frequencies = fromJsonLongArray(json.get("frequencies").getAsJsonArray());
		  // The following reads are not needed as they could well be reconstructed from the info above.
		  //calculatePercentiles();
		  int mode = json.get("mode").getAsInt();
		  modeIndex = Arrays.binarySearch(values, mode);
		  cumulativeProbabilities = fromJsonDoubleArray(json.get("cumulativeProbabilities").getAsJsonArray());
		  percentiles = fromJsonDoubleArray(json.get("percentiles").getAsJsonArray());
		  percentileIndices = fromJsonIntegerArray(json.get("percentileIndices").getAsJsonArray());
		  percentileFrequencyOffsets = fromJsonLongArray(json.get("percentileFrequencyOffsets").getAsJsonArray());
		}
	}

		
	private static double[] fromJsonDoubleArray(Object object) {
		if (object == null) {
			return null;
		}
		else if (object instanceof JsonArray) {
			JsonArray array = (JsonArray) object;
			double[] result = new double[array.size()];
			int nextIndex = 0;
			for (JsonElement e : array) 
				result[nextIndex++] = e.getAsDouble();
			return result;
		}
		else 
			throw new IllegalArgumentException("cannot handle json encoding");
	}
	
	private static int[] fromJsonIntegerArray(Object object) {
		if (object == null) {
			return null;
		}
		else if (object instanceof JsonArray) {
			JsonArray array = (JsonArray) object;
			int[] result = new int[array.size()];
			int nextIndex = 0;
			for (JsonElement e : array) 
				result[nextIndex++] = e.getAsInt();
			return result;
		}
		else 
			throw new IllegalArgumentException("cannot handle json encoding");
	}

	private static long[] fromJsonLongArray(Object object) {
		if (object == null) {
			return null;
		}
		else if (object instanceof JsonArray) {
			JsonArray array = (JsonArray) object;
			long[] result = new long[array.size()];
			int nextIndex = 0;
			for (JsonElement e : array) 
				result[nextIndex++] = e.getAsLong();
			return result;
		}
		else 
			throw new IllegalArgumentException("cannot handle json encoding");
	}
	
	
	
	private static JsonArray toJsonArray(long[] array) {
		JsonArray result = new JsonArray();
		for (long i : array)
			result.add(new JsonPrimitive(i));
		return result;
	}

	private static JsonArray toJsonArray(int[] array) {
		JsonArray result = new JsonArray();
		for (long i : array)
			result.add(new JsonPrimitive(i));
		return result;
	}
	
	
	private static JsonArray toJsonArray(double[] array) {
		JsonArray result = new JsonArray();
		for (double d : array)
			result.add(new JsonPrimitive(d));
		return result;
	}

	@Override
	public long count() {
		return count;
	}	
	
	@Override
	public int mode() {
		if (modeIndex >= 0)
			return values[modeIndex];
		else
			return 0;
	}
	
}
