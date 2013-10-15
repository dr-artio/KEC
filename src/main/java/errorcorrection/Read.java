
package errorcorrection;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;

public class Read extends Sequence {
	public String name;
	ArrayList<Kmer> kmers;
	ArrayList<ErrorRegion> errorRegions;
	int frequency;
	Correction[] corrections;
	int[] NuclType;
	public Read(String s)
	{
		super(s);
		kmers = new ArrayList<Kmer>();
		errorRegions = new ArrayList<ErrorRegion>();
		corrections = new Correction[this.getLength()];
		NuclType = new int[nucl.length()];
	}
	public Read(String s, int f)
	{
		super(s);
		kmers = new ArrayList<Kmer>();
		errorRegions = new ArrayList<ErrorRegion>();
		frequency = f;
		corrections = new Correction[this.getLength()];
		NuclType = new int[nucl.length()];
	}
	public Read(String s, String n)
	{
		super(s);
		name = n;
		kmers = new ArrayList<Kmer>();
		errorRegions = new ArrayList<ErrorRegion>();
		corrections = new Correction[this.getLength()];
		NuclType = new int[nucl.length()];
	}
	public Read(String s, String n, int f)
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
        public double calcEditDistAlign(Read r) throws IOException
                 // <editor-fold defaultstate="collapsed" desc=" DESCRIPTION ">
        {
                FileWriter fw_alin = new FileWriter("allign_input.fas");
                fw_alin.write(">" + this.name + "\n" + this.nucl + "\n");
                fw_alin.write(">" + r.name + "\n" + r.nucl + "\n");
                fw_alin.close();
            
            
                String[] s = new String[] {"Muscle" + File.separator + "muscle", "-in", "allign_input.fas", "-out", "allign_output.fas"};
                    
                ProcessBuilder pb = new ProcessBuilder(s);
                pb.redirectErrorStream(true);
                Process proc = pb.start();
                
                InputStream is = proc.getInputStream();
                InputStreamReader isr = new InputStreamReader(is);
                BufferedReader br = new BufferedReader(isr);

                String ss;
                while ((ss = br.readLine()) != null) {
//                    System.out.println(ss);
                }
                
                if (proc.exitValue()!= 0)
                {
                    System.out.println("Error in allgnment program");
                }
                
 /*               double gapop = 15;
                double gapext = 6.6;
                Runtime run=Runtime.getRuntime();
                String param = "ClustalW2//clustalw2 -INFILE=" +"allign_input.fas" + " -OUTFILE=" + "allign_output.fas";
                param+= " -OUTPUT=FASTA -DNAMATRIX=IUB -GAPOPEN=" + gapop +" -GAPEXT=" + gapext + " -TYPE=DNA -PWDNAMATRIX=IUB -PWGAPOPEN=" + gapop+ " -PWGAPEXT=" + gapext;
                Process p = run.exec(param);
                
                BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
                String s1 = "";
                while ((s1 = stdInput.readLine()) != null)
                {
//                  System.out.println(s1);
                }
                if (p.exitValue()!= 0)
                {
                    System.out.println("Error in allgnment program");
                }*/
            
            
                DataSet alignment = new DataSet("allign_output.fas",'c');
                File f = new File("allign_input.fas");
                f.delete();
                f = new File("allign_output.fas");
                f.delete();
                f = new File("allign_input.dnd");
                f.delete();
                String fir = alignment.reads.get(0).getNucl();
                String sec = alignment.reads.get(1).getNucl();
                int stcomp = 0;
                int endcomp = fir.length() - 1;
                while ((fir.charAt(stcomp) == '-') || (sec.charAt(stcomp) == '-'))
                    stcomp++;
                while ((fir.charAt(endcomp) == '-') || (sec.charAt(endcomp) == '-'))
                    endcomp--;
                int dist = 0;
                for (int i = stcomp; i <= endcomp; i++)
                    if (fir.charAt(i)!=sec.charAt(i))
                        dist++;
                
                double d = ((double) dist)/(endcomp - stcomp + 1);
                return d;
                
        }
        // </editor-fold>
        public int calcEditDistAbsAlign(Read r) throws IOException
                 // <editor-fold defaultstate="collapsed" desc=" DESCRIPTION ">
        {
                FileWriter fw_alin = new FileWriter("allign_input.fas");
                fw_alin.write(">" + this.name + "\n" + this.nucl + "\n");
                fw_alin.write(">" + r.name + "\n" + r.nucl + "\n");
                fw_alin.close();
            
            
                String[] s = new String[] {"Muscle" + File.separator + "muscle", "-in", "allign_input.fas", "-out", "allign_output.fas"};
                    
                ProcessBuilder pb = new ProcessBuilder(s);
                pb.redirectErrorStream(true);
                Process proc = pb.start();
                
                InputStream is = proc.getInputStream();
                InputStreamReader isr = new InputStreamReader(is);
                BufferedReader br = new BufferedReader(isr);

                String ss;
                while ((ss = br.readLine()) != null) {
//                    System.out.println(ss);
                }
                
                if (proc.exitValue()!= 0)
                {
                    System.out.println("Error in allgnment program");
                }
                
 /*               double gapop = 15;
                double gapext = 6.6;
                Runtime run=Runtime.getRuntime();
                String param = "ClustalW2//clustalw2 -INFILE=" +"allign_input.fas" + " -OUTFILE=" + "allign_output.fas";
                param+= " -OUTPUT=FASTA -DNAMATRIX=IUB -GAPOPEN=" + gapop +" -GAPEXT=" + gapext + " -TYPE=DNA -PWDNAMATRIX=IUB -PWGAPOPEN=" + gapop+ " -PWGAPEXT=" + gapext;
                Process p = run.exec(param);
                
                BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
                String s1 = "";
                while ((s1 = stdInput.readLine()) != null)
                {
//                  System.out.println(s1);
                }
                if (p.exitValue()!= 0)
                {
                    System.out.println("Error in allgnment program");
                }*/
            
            
                DataSet alignment = new DataSet("allign_output.fas",'c');
                File f = new File("allign_input.fas");
                f.delete();
                f = new File("allign_output.fas");
                f.delete();
                f = new File("allign_input.dnd");
                f.delete();
                String fir = alignment.reads.get(0).getNucl();
                String sec = alignment.reads.get(1).getNucl();
                int stcomp = 0;
                int endcomp = fir.length() - 1;
                while ((fir.charAt(stcomp) == '-') || (sec.charAt(stcomp) == '-'))
                    stcomp++;
                while ((fir.charAt(endcomp) == '-') || (sec.charAt(endcomp) == '-'))
                    endcomp--;
                int dist = 0;
                for (int i = stcomp; i <= endcomp; i++)
                    if (fir.charAt(i)!=sec.charAt(i))
                        dist++;
                
                return dist;
                
        }
        // </editor-fold>
        public int calcEditDistAbsAlignIgnoreIns(Read r) throws IOException
                 // <editor-fold defaultstate="collapsed" desc=" DESCRIPTION ">
        {
                FileWriter fw_alin = new FileWriter("allign_input.fas");
                fw_alin.write(">" + this.name + "\n" + this.nucl + "\n");
                fw_alin.write(">" + r.name + "\n" + r.nucl + "\n");
                fw_alin.close();
            
            
                String[] s = new String[] {"Muscle" + File.separator + "muscle", "-in", "allign_input.fas", "-out", "allign_output.fas"};
                    
                ProcessBuilder pb = new ProcessBuilder(s);
                pb.redirectErrorStream(true);
                Process proc = pb.start();
                
                InputStream is = proc.getInputStream();
                InputStreamReader isr = new InputStreamReader(is);
                BufferedReader br = new BufferedReader(isr);

                String ss;
                while ((ss = br.readLine()) != null) {
//                    System.out.println(ss);
                }
                
                if (proc.exitValue()!= 0)
                {
                    System.out.println("Error in allgnment program");
                }
                
 /*               double gapop = 15;
                double gapext = 6.6;
                Runtime run=Runtime.getRuntime();
                String param = "ClustalW2//clustalw2 -INFILE=" +"allign_input.fas" + " -OUTFILE=" + "allign_output.fas";
                param+= " -OUTPUT=FASTA -DNAMATRIX=IUB -GAPOPEN=" + gapop +" -GAPEXT=" + gapext + " -TYPE=DNA -PWDNAMATRIX=IUB -PWGAPOPEN=" + gapop+ " -PWGAPEXT=" + gapext;
                Process p = run.exec(param);
                
                BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
                String s1 = "";
                while ((s1 = stdInput.readLine()) != null)
                {
//                  System.out.println(s1);
                }
                if (p.exitValue()!= 0)
                {
                    System.out.println("Error in allgnment program");
                }*/
            
            
                DataSet alignment = new DataSet("allign_output.fas",'c');
                File f = new File("allign_input.fas");
                f.delete();
                f = new File("allign_output.fas");
                f.delete();
                f = new File("allign_input.dnd");
                f.delete();
                String fir = alignment.reads.get(0).getNucl();
                String sec = alignment.reads.get(1).getNucl();
                int stcomp = 0;
                int endcomp = fir.length() - 1;
                while ((fir.charAt(stcomp) == '-') || (sec.charAt(stcomp) == '-'))
                    stcomp++;
                while ((fir.charAt(endcomp) == '-') || (sec.charAt(endcomp) == '-'))
                    endcomp--;
                int dist = 0;
                for (int i = stcomp; i <= endcomp; i++)
                    if ((fir.charAt(i)!=sec.charAt(i)) && (fir.charAt(i) != '-') && (sec.charAt(i) != '-'))
                        dist++;
                
                return dist;
                
        }
        // </editor-fold>
        public double calcEditDistKmer(Read r)
                // <editor-fold defaultstate="collapsed" desc=" DESCRIPTION ">
        {
            HashSet hs1 = new HashSet();
            for (Kmer km : this.kmers)
                hs1.add(km.nucl);
            
            HashSet hs2 = new HashSet();
            for (Kmer km : r.kmers)
                hs2.add(km.nucl);
            
            int n1 = hs1.size();
            int n2 = hs2.size();
            
            hs1.retainAll(hs2);
            int m = hs1.size();
            double d = 1-((double )m)/Math.min(n1, n2);
            return d;
        }
        // </editor-fold>
        public int getFreq()
        {
            return frequency;
        }
        public double calcEditDistAlign(Read r, int gapop, int gapext) throws IOException
                 // <editor-fold defaultstate="collapsed" desc=" DESCRIPTION ">
        {

		DNASequence target = new DNASequence(this.nucl,
				AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence query = new DNASequence(r.nucl,
				AmbiguityDNACompoundSet.getDNACompoundSet());
 
		SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
 
		SimpleGapPenalty gapP = new SimpleGapPenalty();
		gapP.setOpenPenalty((short)gapop);
		gapP.setExtensionPenalty((short)gapext);
 
		SequencePair<DNASequence, NucleotideCompound> psa =
				Alignments.getPairwiseAlignment(query, target,
						PairwiseSequenceAlignerType.GLOBAL, gapP, matrix);
                
                
                int st = 1;
                while (psa.hasGap(st))
                    st++;
                int end = psa.getLength();
                while (psa.hasGap(end))
                    end--;
                
                int d = 0;
                int m = 0;
                for (int i = st; i <= end; i++)
                {
//                    if ((psa.getCompoundAt(1, i) != psa.getCompoundAt(2, i)) && !(psa.hasGap(i)))
                    if (psa.getCompoundAt(1, i) != psa.getCompoundAt(2, i))
                    {
//                        System.out.println(i);
                        d++;
                    }
                    if (!psa.hasGap(i))
                        m++;
                }
                
/*                System.out.println(psa);
                System.out.println("Start: " + st);
                System.out.println("End: " + end);
                System.out.println("nIdent: " + psa.getNumIdenticals());
                System.err.println("d: " + d);
                System.out.println(((double) d)/(end-st+1));*/
               
                
                
                return ((double) d)/(end-st+1);
//                 return ((double) d)/m;
//                 return d;
                
        }
        // </editor-fold>
        public double calcEditDistAbsAlign(Read r, int gapop, int gapext) throws IOException
                 // <editor-fold defaultstate="collapsed" desc=" DESCRIPTION ">
        {

		DNASequence target = new DNASequence(this.nucl,
				AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence query = new DNASequence(r.nucl,
				AmbiguityDNACompoundSet.getDNACompoundSet());
 
		SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
 
		SimpleGapPenalty gapP = new SimpleGapPenalty();
		gapP.setOpenPenalty((short)gapop);
		gapP.setExtensionPenalty((short)gapext);
 
		SequencePair<DNASequence, NucleotideCompound> psa =
				Alignments.getPairwiseAlignment(query, target,
						PairwiseSequenceAlignerType.GLOBAL, gapP, matrix);
                
                
                int st = 1;
                while (psa.hasGap(st))
                    st++;
                int end = psa.getLength();
                while (psa.hasGap(end))
                    end--;
                
                int d = 0;
                for (int i = st; i <= end; i++)
                    if ((psa.getCompoundAt(1, i) != psa.getCompoundAt(2, i)) && !(psa.hasGap(i)))
                    {
//                        System.out.println(i);
                        d++;
                    }
                
/*                System.out.println(psa);
                System.out.println(st);
                System.out.println(end);
                System.out.println(psa.getNumIdenticals());
                System.out.println(((double) d)/(end-st+1));*/
               
                
                
                return ((double) d)/(end-st+1);
                
        }
}