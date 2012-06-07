package net.malariagen.gatk.genotyper.models;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;
import java.lang.annotation.ElementType;

@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.TYPE)
public @interface GModel {
	public final static String NO_NAME = "<no-name>";
	String fullName() default NO_NAME;
	String shortName() default NO_NAME;
	String description() default "";
	String[] aliases() default {};
}
