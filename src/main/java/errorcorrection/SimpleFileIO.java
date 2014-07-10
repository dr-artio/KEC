package ErrorCorrection;


import static java.lang.String.format;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Type;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.StringTokenizer;

/**
 * Class provide an interface for 
 * reading/writing collections and files
 * 
 * @author ava
 */
public class SimpleFileIO {
	/**
	 * Save sorted map using custom comparator for sorting
	 * 
	 * @param filename
	 *            Location and file name for saving
	 * @param map
	 *            for saving
	 * @param format
	 *            for toString()
	 * @param args
	 *            additional arguments
	 */
	public static <TKey, TValue> void saveMap(final String filename,
			final Map<TKey, TValue> map, final String format,
			final Object... args) {
		save(filename, getEntryList(map, format), args);
	}

	/**
	 * Save sorted map using custom comparator for sorting
	 * 
	 * @param filename
	 *            Location and file name for saving
	 * @param map
	 *            for saving
	 * @param comparator
	 *            for sorting
	 * @param format
	 *            for toString()
	 * @param args
	 *            additional arguments
	 */
	public static <TKey, TValue> void saveMap(final String filename,
			final Map<TKey, TValue> map,
			Comparator<Entry<TKey, TValue>> comparator, final String format,
			Object... args) {
		List<Entry<TKey, TValue>> tmp = getEntryList(map, 
				format == null || format.equals("") ? "%s:%s%n" : format);
		Collections.sort(tmp, comparator);
		save(filename, tmp, args);
	}

	/**
	 * Save collection into file. Carefully toString() method should be
	 * correctly overridden
	 * 
	 * @param filename
	 * @param objs
	 * @param args
	 */
	public static <T> void save(final String filename,
			final Collection<T> objs, Object... args) {
		try {
			save(new FileOutputStream(format(filename, args)),objs, args);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Save collection into stream. Carefully toString() method should be
	 * correctly overridden
	 * 
	 * @param stream
	 * @param objs
	 * @param args
	 */
	public static <T> void save(final OutputStream stream,
			final Collection<T> objs, Object... args) throws IOException{
		PrintWriter out = null;
		try {
			out = new PrintWriter(stream);
			for (final T o : objs)
				out.write(o.toString()+"\n");
		} finally {
			if (out != null)
				out.close();
		}
	}

	/**
	 * Read file as a collection of strings
	 * 
	 * @param filename
	 * @return
	 * @throws IOException
	 */
	public static Collection<String> readStrings(final String filename)
			throws IOException {
		Collection<String> result = new ArrayList<String>();
		String line;
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(filename));
			while ((line = reader.readLine()) != null) {
				result.add(line.trim());
			}
			return result;
		} finally {
			if (reader != null)
				reader.close();
		}
	}
        
        /**
	 * Read file as a collection of strings
	 * 
	 * @param filename
	 * @return
	 * @throws IOException
	 */
	public static<T1,T2> Map<T1,T2> read(final String filename)
			throws IOException, NoSuchMethodException, InstantiationException, IllegalAccessException, IllegalArgumentException, InvocationTargetException {
            Map<T1, T2> map = new HashMap<T1, T2>();
            Collection<String> strs = readStrings(filename);
            Type[] genericTypes = SimpleFileIO.class.getDeclaredMethod("read", String.class).getGenericParameterTypes();
            Constructor<? extends Type> keyConstructor = genericTypes[0].getClass().getConstructor(String.class);
            Constructor<? extends Type> valueConstructor = genericTypes[1].getClass().getConstructor(String.class);
            for(String s : strs){
                StringTokenizer st = new StringTokenizer(s," ");
                if (st.countTokens()>1){
                    map.put((T1)keyConstructor.newInstance(st.nextToken()),(T2)valueConstructor.newInstance(st.nextToken()));
                }   
            }
            return map;
        }
        
        
	/**
	 * @param map
	 *            for saving
	 * @param format
	 *            string for toString method been overridden
	 *            first parameter is TKey, second - TValue
	 * @return List of Entries
	 */
	@SuppressWarnings("serial")
	private static <TKey, TValue> List<Entry<TKey, TValue>> getEntryList(
			final Map<TKey, TValue> map, final String format) {
		if (map == null) return new ArrayList<Entry<TKey, TValue>>();
		List<Entry<TKey, TValue>> tmp = new ArrayList<Entry<TKey, TValue>>(
				map.size());
		for (Entry<TKey, TValue> k : map.entrySet()) {
			tmp.add(new AbstractMap.SimpleEntry<TKey, TValue>(k) {
				@Override
				public String toString() {
					return format(format, getKey(), getValue());
				}
			});
		}
		return tmp;
	}

}
