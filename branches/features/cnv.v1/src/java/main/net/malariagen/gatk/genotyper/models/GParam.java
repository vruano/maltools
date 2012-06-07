package net.malariagen.gatk.genotyper.models;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.FIELD)
public @interface GParam {
	String shortName() default NO_NAME;
	String fullName() default NO_NAME;
	int offset() default -1;
	boolean required() default false;
	String description() default "";
	
	public static String NO_NAME = "<no-name>";
	
}
