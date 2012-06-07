package net.malariagen.utils;

import java.util.LinkedList;
import java.util.List;


/**
 * Non-thread safe but fast implementation for an object Pool.
 * 
 * @author valentin
 *
 * @param <T>
 */
public class SequentialPool<T> implements Pool<T> {

	private final Factory<T> factory;
	private final List<T> queue;
	
	public SequentialPool() {
		this(null);
	}
	
	public SequentialPool(Factory<T> factory) {
		this.factory = factory;
		queue = new LinkedList<T>();
	}
	
	@Override
	public T borrow() {
		if (queue.isEmpty()) {
			if (factory != null) return factory.newInstance();
			else return null;
		}
		else 
			return queue.remove(0);
	}
	
	@Override
	public void restore(T o) {
		queue.add(o);
	}
	
}
