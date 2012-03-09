package net.malariagen.utils;

public interface Pool<T> {

	T borrow();

	void restore(T o);

}
