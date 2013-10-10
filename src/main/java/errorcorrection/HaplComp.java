package errorcorrection;

import java.util.Comparator;

public class HaplComp implements Comparator<Haplotype>{
		public int compare(Haplotype r1, Haplotype r2) {
			return -(r1.frequency - r2.frequency);
		}
}
