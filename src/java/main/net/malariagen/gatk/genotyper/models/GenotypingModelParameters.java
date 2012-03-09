package net.malariagen.gatk.genotyper.models;

import java.io.File;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.sting.utils.collections.Pair;

public class GenotypingModelParameters<M extends GenotypingModel> {

	private Class<M> modelClass;

	private List<Field> paramFields;
	private Map<String, Field> fieldByName;

	public GenotypingModelParameters(Class<M> modelClass) {
		if (modelClass == null)
			throw new IllegalArgumentException();
		this.modelClass = modelClass;
		paramFields = new ArrayList<Field>(10);
		fieldByName = new HashMap<String, Field>(10);
		findParametersNamePass();
		findParametersOffsetPass();
		findParameterRemainderPass(modelClass, 0);
	}

	public Class<M> getModelClass() {
		return modelClass;
	}

	private void findParametersNamePass() {
		Class<? super M> cls = modelClass;
		while (cls != null) {
			Field[] fields = cls.getDeclaredFields();
			for (Field f : fields) {
				GParam anno = f.getAnnotation(GParam.class);
				if (anno == null)
					continue;
				f.setAccessible(true);
				if (anno.shortName() != null) {
					if (!fieldByName.containsKey(anno.shortName()))
						fieldByName.put(anno.shortName(), f);
				}
				if (anno.fullName() != null) {
					if (!fieldByName.containsKey(anno.fullName()))
						fieldByName.put(anno.fullName(), f);
				}
			}
			cls = cls.getSuperclass();
		}
	}

	private void findParametersOffsetPass() {
		Class<? super M> cls = modelClass;
		while (cls != null) {
			Field[] fields = cls.getDeclaredFields();
			for (Field f : fields) {
				GParam anno = f.getAnnotation(GParam.class);
				if (anno == null)
					continue;
				if (anno.offset() >= 0) {
					int offset = anno.offset();
					if (paramFields.size() <= offset) {
						for (int i = paramFields.size(); i < offset; i++)
							paramFields.add(null);
						paramFields.add(f);
					} else if (paramFields.get(offset) == null)
						paramFields.set(offset, f);
				}
			}
			cls = cls.getSuperclass();
		}
	}

	private int findParameterRemainderPass(Class<? super M> cls, int fromIndex) {
		if (cls.getSuperclass() != null)
			fromIndex = findParameterRemainderPass(cls.getSuperclass(),
					fromIndex);
		Field[] fields = cls.getDeclaredFields();
		for (Field f : fields) {
			GParam anno = f.getAnnotation(GParam.class);
			if (anno == null)
				continue;
			if (anno.offset() >= 0)
				continue;
			if (!anno.shortName().equals(GParam.NO_NAME)) {
				if (anno.fullName() != null) {
					if (!f.equals(fieldByName.get(anno.shortName()))
							&& !f.equals(fieldByName.get(anno.fullName())))
						continue;
				} else if (!f.equals(fieldByName.get(anno.shortName())))
					continue;
			} else if (anno.fullName() != null
					&& !f.equals(fieldByName.get(anno.fullName())))
				continue;
			while (paramFields.size() > fromIndex
					&& paramFields.get(fromIndex) != null)
				fromIndex++;
			if (paramFields.size() == fromIndex)
				paramFields.add(f);
			else
				paramFields.set(fromIndex, f);
			fromIndex++;
		}
		return fromIndex;

	}

	public String getModelDisplayName() {
		GModel anno = modelClass.getAnnotation(GModel.class);
		if (!GModel.NO_NAME.equals(anno.shortName()))
			return anno.shortName();
		if (!GModel.NO_NAME.equals(anno.fullName()))
			return anno.fullName();
		if (anno.aliases().length > 0)
			return anno.aliases()[0];
		return modelClass.getSimpleName();
	}

	public void applyParameters(M model, List<Pair<String,String>> params) {
		if (!model.getClass().equals(modelClass))
			throw new IllegalArgumentException("you must provide a model of the same concrete class as this parameter set's");
		int nextOffset = 0;
		Set<Field> provided = new HashSet<Field>();
		for (Pair<String, String> p : params) {
			if (p.first != null) {
				Field f = fieldByName.get(p.first);
				if (f == null)
					throw new GenotypingModelException("there is not such paramater with name '" + p.first + "' in model '" + getModelDisplayName() + "'");
				int fOffset = f.getAnnotation(GParam.class).offset();
				if (fOffset == nextOffset)
					nextOffset++;
				if (provided.contains(f)) 
					throw new GenotypingModelException("parameter '" + displayName(f) + "' provided more than once in model '" + getModelDisplayName() + "'");	
				applyParameter(model,f,p.second);
				provided.add(f);
			}
			else {
				Field f;
				while (paramFields.size() > nextOffset) {
					f = paramFields.get(nextOffset);
					if (!provided.contains(f)) {
						applyParameter(model,f,p.second);
						provided.add(f);
						break;
					}
					nextOffset++;
				}
				if (paramFields.size() <= nextOffset) 
					throw new GenotypingModelException("too many parameters '" + params.size() + " < " + paramFields.size() + "' specified for model '" + getModelDisplayName() + "'");
			}
		}
	}

	private String displayName(Field f) {
		GParam anno = f.getAnnotation(GParam.class);
		if (!anno.shortName().equals(GParam.NO_NAME))
			return anno.shortName();
		else if (!anno.fullName().equals(GParam.NO_NAME))
			return anno.fullName();
		else if (anno.offset() >= 0) {
			int offset = anno.offset();
			int pos = offset + 1;
			if (pos == 1)
				return "1st";
			else if (pos == 2)
				return "2nd";
			else if (pos == 3)
				return "3rd";
			else if (pos % 100 > 20) {
				if (pos % 10 == 1)
					return pos + "st";
				else if (pos % 10 == 2)
					return pos + "nd";
				else if (pos % 10 == 3)
					return pos + "rd";
				else
					return pos + "th";
			} else
				return pos + "th";
		} else
			return f.getName();
	}
	
	private void applyParameter(M model, Field f, String value) {
		try {
			Class<?> type = f.getType();
			if (type.isPrimitive()) {
				if (Integer.TYPE.equals(type)) {
					f.setInt(model, Integer.parseInt(value));
				} else if (Boolean.TYPE.equals(type)) {
					f.setBoolean(model, Boolean.parseBoolean(value));
				} else if (Long.TYPE.equals(type)) {
					f.setLong(model, Long.parseLong(value));
				} else if (Double.TYPE.equals(type)) {
					f.setDouble(model, Double.parseDouble(value));
				} else if (Float.TYPE.equals(type)) {
					f.setFloat(model, Float.parseFloat(value));
				} else if (Short.TYPE.equals(type)) {
					f.setShort(model, Short.parseShort(value));
				} else if (Byte.TYPE.equals(type)) {
					f.setByte(model, Byte.parseByte(value));
				} else if (Character.TYPE.equals(type)) {
					if (value.length() != 1)
						throw new GenotypingModelException(
								"char-type model parameter '"
										+ displayName(f)
										+ "' in model '"
										+ getModelDisplayName()
										+ "' only accept single character values");
					f.setChar(model, value.charAt(0));
				}
			} else if (type.getName().equals(String.class.getName())) {
				f.set(model, value);
			} else if (type.getName().equals(File.class.getName())) {
				f.set(model, new File(value));
			} else if (type.isEnum()) {
				@SuppressWarnings("unchecked")
				Class<? extends Enum<?>> enumCls = (Class<? extends Enum<?>>) type;
				Object v;
				try {
					v = enumCls.getMethod("valueOf", String.class)
							.invoke(value);
				} catch (IllegalArgumentException e) {
					throw new GenotypingModelException(
							"invalid model parameter type '"
									+ type.getName()
									+ "' for parameter '"
									+ displayName(f)
									+ "' in model '"
									+ getModelDisplayName() + "'", e);
				} catch (Exception e) {
					throw new GenotypingModelException(
							"unexpected Enum.valueOf failure at model parameter '"
									+ displayName(f)
									+ "' in model '"
									+ getModelDisplayName() + "'", e);
				}
				f.set(f, v);
			} else if (Number.class.isAssignableFrom(type)) {
				if (Double.class.isAssignableFrom(type)) {
					f.set(model,
							new Double(Double.parseDouble(value)));
				} else if (Integer.class.isAssignableFrom(type)) {
					f
							.set(model, new Integer(Integer.parseInt(value)));
				} else if (Short.class.isAssignableFrom(type)) {
					f.set(model, new Short(Short.parseShort(value)));
				} else if (Byte.class.isAssignableFrom(type)) {
					f.set(model, new Byte(Byte.parseByte(value)));
				} else if (Float.class.isAssignableFrom(type)) {
					f.set(model, new Float(Float.parseFloat(value)));
				}
			} else if (Boolean.class.isAssignableFrom(type)) {
				f.set(model, Boolean.valueOf(value));
			} else if (Character.class.isAssignableFrom(type)) {
				if (value.length() != 1)
					throw new GenotypingModelException(
							"char-type model parameter '"
									+ displayName(f)
									+ "' in model '"
									+ getModelDisplayName() + "' only accept single character values");
				f.set(model, Character.valueOf(value.charAt(0)));
			} else {
				throw new GenotypingModelException(
						"invalid model parameter type '"
								+ type.getName()
								+ "' for parameter '"
								+ displayName(f)
								+ "' in model '"
								+ getModelDisplayName() + "'");
			}

		} catch (NumberFormatException e) {
			throw new GenotypingModelException(
					"invalid parameter value for parameter '"
							+ displayName(f)
							+ "' in model '"
							+ getModelDisplayName(), e);
		} catch (IllegalArgumentException e) {
			throw new GenotypingModelException(
					"invalid parameter value for parameter '"
							+ displayName(f)
							+ "' in model '"
							+ getModelDisplayName(), e);
		} catch (IllegalAccessException e) {
			throw new GenotypingModelException(
					"invalid parameter definition for parameter '"
							+ displayName(f)
							+ "' in model '"
							+ getModelDisplayName(), e);
		}

	}

}