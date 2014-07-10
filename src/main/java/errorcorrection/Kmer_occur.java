package ErrorCorrection;

public class Kmer_occur {
	Read r;
	int readNumber;
	int stPos;
	int length;
	Kmer_occur(Read r1, int rn, int i, int k)
	{
		r = r1;
		readNumber = rn;
		stPos = i;
		length = k;
	}
	void println()
	{
		System.out.println(r.name + " " + stPos + " " + length);
	}

}
