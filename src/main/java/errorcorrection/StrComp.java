package ErrorCorrection;

import java.util.Comparator;

public class StrComp implements Comparator<Read> {

	@Override
	public int compare(Read r1, Read r2) {
		return -(r1.nucl.length() - r2.nucl.length());
	}

}
