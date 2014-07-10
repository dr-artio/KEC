package ErrorCorrection;

import java.io.IOException;

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
		if (nucl.contains(r.nucl))
			return true;
		else
			return false;
	}
        boolean containsReadAlign(Read r, int gapop, int gapext) throws IOException
	{
		Read s = new Read(this.nucl);
                if (s.calcEditDistAbsAlign(r, gapop, gapext) == 0)
                    return true;
                else
                    return false;
	}
}
