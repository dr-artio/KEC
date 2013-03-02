package errorcorrection;

public class Haplotype extends Sequence {
	int frequency;
	String name;
	double concentration;
	String clone;
        
	Haplotype(Read r)
	{
		nucl = r.nucl;
		name = r.name;
		frequency = r.frequency;
		concentration = r.frequency;
		clone = "unknown";
	}
        
	boolean containsRead(Read r)
	{
		return nucl.contains(r.nucl);
	}
}
