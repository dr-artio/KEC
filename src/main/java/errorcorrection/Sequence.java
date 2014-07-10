package ErrorCorrection;
import java.util.*;

public class Sequence {
	String nucl;
	public Sequence()
	{
		nucl = new String();
	}
	public Sequence(int l)
	{
		nucl = new String();
	}
	public Sequence(String s)
	{
		nucl = s;
	}
	int getLength()
	{
		return nucl.length();
	}
	char getNucl(int i)
	{
		return nucl.charAt(i);
	}
	public String getNucl(int i, int j)
	{
		return nucl.substring(i, j+1);
	}
	public String getNucl()
	{
		return nucl;
	}
	void removeNucl(int i, int j)
	{
		StringBuffer s = new StringBuffer(nucl);
		s.delete(i, j+1);
		nucl = new String(s);
	}
        void removeNucl(int i)
	{
		StringBuffer s = new StringBuffer(nucl);
                s.deleteCharAt(i);
		nucl = new String(s);
	}
	void replaceNucl(int i, char c)
	{
		StringBuffer s = new StringBuffer(nucl);
		s.setCharAt(i, c);
		nucl = new String(s);
	}
	
	
}
