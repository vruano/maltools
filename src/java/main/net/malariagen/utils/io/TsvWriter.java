package net.malariagen.utils.io;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

public class TsvWriter {

	private Writer writer;
	private String sep;
	private StringBuffer sb;
	
	public TsvWriter(Writer w, String sep) {
		if (w == null)
			throw new IllegalArgumentException("the provided writer cannot be null");
		if (sep == null)
			throw new IllegalArgumentException("the provided separator cannot be null");
		if (sep.indexOf("\n") != -1)
			throw new IllegalArgumentException("the separator cannot contain the new line character");
		this.writer = w;
		this.sep = sep;
		this.sb = new StringBuffer(100);
	}
	
	public TsvWriter(Writer w) {
		this(w,"\t");
	}
	
	public void writeLine(Object ... elements) throws IOException {
		sb.setLength(0);
		for (Object e : elements) {
			if (e == null)
				throw new IllegalArgumentException("no null elements allowed");
			String s = e.toString();
			if (s.indexOf(sep) != -1)
				throw new IllegalArgumentException("elements cannot contain the separator");
			if (s.indexOf("\n") != -1)
				throw new IllegalArgumentException("elements cannot contain the new line character");
			sb.append(s).append(sep);
		}
		if (sb.length() > 0)
			sb.setLength(sb.length() - sep.length());
		sb.append('\n');
		writer.write(sb.toString());
	}
	
	public void writeLine(int ... values) throws IOException {
		sb.setLength(0);
		for (int v : values) {
			if (sep.indexOf("" + v) != -1)
				throw new IllegalArgumentException("elements cannot contain the separator");
			sb.append(v).append(sep);
		}
		if (sb.length() > 0)
			sb.setLength(sb.length() - sep.length());
		sb.append('\n');
		writer.write(sb.toString());
	}

	public void writeLine(long ... values) throws IOException {
		sb.setLength(0);
		for (long v : values) {
			if (sep.indexOf("" + v) != -1)
				throw new IllegalArgumentException("elements cannot contain the separator");
			sb.append(v).append(sep);
		}
		if (sb.length() > 0)
			sb.setLength(sb.length() - sep.length());
		sb.append('\n');
		writer.write(sb.toString());
	}
	
	
	public void writeComment(String comment) throws IOException {
		if (comment == null)
			throw new IllegalArgumentException("comment cannot be null");
		int firstNl = comment.indexOf('\n');
		if (firstNl == -1 || firstNl == comment.length() - 1) {
			if (!comment.startsWith("#"))
				writer.write("# ");
			writer.write(comment);
			if (firstNl == -1)
				writer.write("\n");
		}
		else 
			for (String s : comment.split("\\n"))
				writeComment(s);
	}
	
	public void close() throws IOException {
		writer.close();
	}
	
	public void flush() throws IOException {
		writer.flush();
	}
	
	
	

}
