package net.malariagen.gatk.genotyper.models;

import java.util.ArrayList;
import java.util.List;

import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.reflections.Reflections;

public class GenotypingModelUtils {

	private static final Reflections REFLECTIONS = new Reflections(
			"net.malariagen");

	public synchronized static <M extends GenotypingModel> M getModelInstance(String spec) {
		StringBuffer nameBuffer = new StringBuffer(100);
		List<Pair<String, String>> params = new ArrayList<Pair<String, String>>(
				10);
		parseModelSpec(spec, nameBuffer, params);
		String name = nameBuffer.toString();
		Class<M> modelClass = findModelClassBySpecName(name);
		if (modelClass == null)
			throw new UserException("there is no model with name '" + name
					+ "' in specification '" + spec
					+ "'. Perhaps you meant one of the following: "
					+ allModelNameString());
		GenotypingModelParameters<M> gmp = new GenotypingModelParameters<M>(modelClass);
		M result;
		try {
			result = modelClass.newInstance();
		} catch (InstantiationException e) {
			throw new GenotypingModelException(
					"could not instanciate the genotyping-model class "
							+ modelClass.getName(), e);
		} catch (IllegalAccessException e) {
			throw new GenotypingModelException(
					"could not instanciate the genotyping-model class "
							+ modelClass.getName(), e);
		}
		gmp.applyParameters(result, params);
		return result;
	}


	private static String allModelNameString() {
		List<String> names = new ArrayList<String>(100);
		for (Class<?> cls : REFLECTIONS.getTypesAnnotatedWith(GModel.class)) {
			if (!GenotypingModel.class.isAssignableFrom(cls))
				continue;
			GModel anno = cls.getAnnotation(GModel.class);
			String bestName = null;
			if (anno.shortName() != null) {
				bestName = anno.shortName();
				break;
			} else if (anno.fullName() != null) {
				bestName = anno.fullName();
				break;
			} else if (anno.aliases() != null && anno.aliases().length > 0) {
				bestName = anno.aliases()[0];
			}
			if (bestName != null)
				names.add(bestName);
		}
		if (names.size() == 0) {
			return "<NONE-AVIALABLE>";
		} else if (names.size() == 1) {
			return "'" + names.get(0) + "'";
		} else {
			StringBuffer result = new StringBuffer(names.size() * 100);
			for (int i = 0; i < names.size() - 2; i++) {
				result.append('\'').append(names.get(i)).append("\', ");
			}
			result.append('\'').append(names.get(names.size() - 2))
					.append("\' and '").append(names.get(names.size() - 1))
					.append('\'');
			return result.toString();
		}
	}

	private static <M extends GenotypingModel> Class<M> findModelClassBySpecName(
			String name) {
		for (Class<?> cls : REFLECTIONS.getTypesAnnotatedWith(GModel.class)) {
			if (!GenotypingModel.class.isAssignableFrom(cls))
				continue;
			@SuppressWarnings("unchecked")
			Class<M> gmCls = (Class<M>) cls;
			GModel anno = cls.getAnnotation(GModel.class);
			if (name.equals(anno.shortName())) {
				return gmCls;
			} else if (name.equals(anno.fullName())) {
				return gmCls;
			} else if (anno.aliases() != null) {
				for (String alias : anno.aliases()) {
					if (name.equals(alias)) {
						return gmCls;
					}
				}
			}
		}
		return null;
	}

	private static void parseModelSpec(String spec, StringBuffer nameBuffer,
			List<Pair<String, String>> params) {
		if (spec.matches("^\\s*[a-zA-Z_][a-zA-Z_0-9]*\\s*$")) {
			nameBuffer.append(spec);
		}
		else if (spec.matches("^\\s*[a-zA-Z_][a-zA-Z_0-9]*\\s*\\[[^\\]]*\\]\\s*$")) {
			parseParametricModelSpec(spec,'[',']',nameBuffer,params);
		}
		else if (spec.matches("^\\s*[a-zA-Z_][a-zA-Z_0-9]*\\s*\\([^\\(]*\\)\\s*$")) {
			parseParametricModelSpec(spec,'(',')',nameBuffer,params);
		}
		else if (spec.matches("^\\s*[a-zA-Z_][a-zA-Z_0-9]*\\s*\\{[^\\}]*\\}\\s*$")) {
			parseParametricModelSpec(spec,'{','}',nameBuffer,params);
		}
		else if (spec.matches("^\\s*[a-zA-Z_][a-zA-Z_0-9]*\\s*\\/[^\\/]*\\/\\s*$")) {
			parseParametricModelSpec(spec,'/','/',nameBuffer,params);
		}
		else {
			throw new GenotypingModelException("Illegal genotyping model specification '" + spec + "'");
		}
	}
	
	private static void parseParametricModelSpec(String spec,  char openParams, char closeParams, StringBuffer nameBuffer,List<Pair<String,String>> params) {
		int openParamsPos = spec.indexOf(openParams);
		int closeParamsPos = spec.lastIndexOf(closeParams);
		nameBuffer.append(spec.substring(0,openParamsPos).trim());
		String[] paramStrs = spec.substring(openParamsPos + 1,closeParamsPos).split(",");
		if (paramStrs.length == 1 && paramStrs[0].trim().isEmpty())
			return;
		for (String ps : paramStrs) {
			int equalIndex = ps.indexOf("=");
			if (equalIndex == -1) {
				params.add(new Pair<String,String>(null,ps.trim()));
			}
			else {
				params.add(new Pair<String,String>(ps.substring(0,equalIndex).trim(),ps.substring(equalIndex + 1)));
			}
		}
	}

}
