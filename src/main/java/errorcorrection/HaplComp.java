package ErrorCorrection;

import java.util.Comparator;

public class HaplComp implements Comparator<Haplotype>{
		public int compare(Haplotype r1, Haplotype r2) {
                    Double d1 = r1.concentration;
                    Double d2 = r2.concentration;
                    return -d1.compareTo(d2);
		}
}
