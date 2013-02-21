package ErrorCorrection;
import java.util.*;
import java.math.*;

public class Kmer extends Sequence{
	int kCount;
//	int kCount1;
	int StartPos;
//	LinkedList<Integer> reads;
	double clusterParam1;
//	double clusterParam2;
	int getCount()
	{
		return kCount;
	}
	Kmer(String s, int st)
	{
		super(s);
		kCount = 1;
		StartPos = st;
		clusterParam1 = -1;
	}
	void setKCount(int i)
	{
		kCount = i;
	}
	int getKCount()
	{
		return kCount;
	}
	int getStartPos()
	{
		return StartPos;
	}
	void setClusterParam (double c1)
	{
		clusterParam1 = c1;
//		clusterParam2 = c2;
	}
	public boolean equals(Object ob)
	{
		return nucl.equals(((Kmer) ob).nucl);
	}
}
