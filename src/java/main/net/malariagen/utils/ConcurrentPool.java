package net.malariagen.utils;

import java.util.concurrent.ConcurrentLinkedQueue;

/**
 * Thread-safe but slower implementation for an object pool.
 * 
 * @author valentin
 *
 * @param <T>
 */
public class ConcurrentPool<T> implements Pool<T> {

	private final Factory<T> factory;
	private final ConcurrentLinkedQueue<T> queue;
	
	
	public ConcurrentPool() {
		this(null);
	}
	
	public ConcurrentPool(Factory<T> factory) {
		this.factory = factory;
		queue = new ConcurrentLinkedQueue<T>();
	}
	
	@Override
	public T borrow() {
		T result = queue.poll();
		if (result == null) {
			if (factory != null) return factory.newInstance();
			else return null;
		}
		else 
			return result;
	}
	
	@Override
	public void restore(T o) {
		queue.add(o);
	}
}
