
package errorcorrection;

import java.util.*;

public class Read extends Sequence {
	public String name;
	ArrayList<Kmer> kmers;
	ArrayList<ErrorRegion> errorRegions;
	int frequency;
	Correction[] corrections;
	int[] NuclType;
	Read(String s)
	{
		super(s);
		kmers = new ArrayList<Kmer>();
		errorRegions = new ArrayList<ErrorRegion>();
		corrections = new Correction[this.getLength()];
		NuclType = new int[nucl.length()];
	}
	Read(String s, int f)
	{
		super(s);
		kmers = new ArrayList<Kmer>();
		errorRegions = new ArrayList<ErrorRegion>();
		frequency = f;
		corrections = new Correction[this.getLength()];
		NuclType = new int[nucl.length()];
	}
	Read(String s, String n)
	{
		super(s);
		name = n;
		kmers = new ArrayList<Kmer>();
		errorRegions = new ArrayList<ErrorRegion>();
		corrections = new Correction[this.getLength()];
		NuclType = new int[nucl.length()];
	}
	Read(String s, String n, int f)
	{
		super(s);
		name = n;
		frequency = f;
		kmers = new ArrayList<Kmer>();
		errorRegions = new ArrayList<ErrorRegion>();
		corrections = new Correction[this.getLength()];
		NuclType = new int[nucl.length()];
	}
	void addKmer(Kmer k)
	{
		kmers.add(k);
	}
	public int getLength()
	{
		return nucl.length();
	}
	void addErrorRegion(ErrorRegion e)
	{
		errorRegions.add(e);
	}
        Read RevComp()
        {
            String rev_ref = "";
            for (int i = this.getLength()-1; i >= 0; i--)
            {
                if (this.nucl.charAt(i) == 'A')
                    rev_ref += 'T';
                if (this.nucl.charAt(i) == 'T')
                    rev_ref += 'A';
                if (this.nucl.charAt(i) == 'G')
                    rev_ref += 'C';
                if (this.nucl.charAt(i) == 'C')
                    rev_ref += 'G';
            }
            return new Read(rev_ref);
        }
        Read()
                // <editor-fold defaultstate="collapsed" desc=" DESCRIPTION ">
        {
            
        }
        // </editor-fold>
}