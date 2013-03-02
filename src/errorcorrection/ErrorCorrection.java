package errorcorrection;

import java.io.File;
import java.io.IOException;
import net.sourceforge.argparse4j.*;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.*;

public class ErrorCorrection {

    private static final String K_PARAMETER = "K";
    private static final String I_PARAMETER = "I";
    private static final String L_PARAMETER = "L";
    private static final String ALIGN_PARAMETER = "align";
    private static final String READS_PARAMETER = "reads";
    private static final String DOMINGEN_PARAMETER = "domgen";
    private static final String DOMINPOSTPROC_PARAMETER = "dompostproc";
    static File env_path;

    static {
        env_path = new File(".");
    }

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) 
            throws IOException, InterruptedException {
        ArgumentParser parser = ArgumentParsers.newArgumentParser("KEC")
                .description("Correct genomic sequence data using k-mer count information.");

        int dominparamgen = 30;
        int dominparampostpr = 30;

        int k = 25;
        int nIter = 3;
        boolean toFindHapl = false;
        int errorsseglen = 0;
        File fl = null;

        parser.addArgument(READS_PARAMETER)
                .metavar("ReadsFile")
                .help("File containing raw sequencing data"
                + " file with extension .fasta (.fas) "
                + "reads in fasta format")
                .type(File.class);

        parser.addArgument("-k").dest(K_PARAMETER)
                .metavar("K")
                .setDefault(k)
                .type(Integer.class)
                .help("Parameter k - the size of substrings (k-mers) "
                + "Default: " + k + ")");

        parser.addArgument("-i").dest(I_PARAMETER)
                .metavar("Iterations")
                .setDefault(nIter)
                .type(Integer.class)
                .help("Number of iterations of cleaning"
                + " in error correction procedure. If it is "
                + "greater than 3 than influence is not significant "
                + "(Default: " + nIter + ")");

        parser.addArgument(new String[]{"-a", "--align"}).dest(ALIGN_PARAMETER)
                .action(Arguments.storeTrue())
                .setDefault(false)
                .type(Boolean.class)
                .help("Enable using of CLustalW for multiple sequence "
                + "alignment for additional correction procedure (Default: "
                + "do not align)");

        parser.addArgument("-l").dest(L_PARAMETER)
                .metavar("Num of zeros")
                .setDefault(errorsseglen)
                .type(Integer.class)
                .help("Number of consecutef zeros for finding threshold. "
                + "Using Poisson distribution if 0 (Default: " + errorsseglen + ")");

        parser.addArgument("-dg").dest(DOMINGEN_PARAMETER)
                .metavar("dominparamgen")
                .setDefault(dominparamgen)
                .type(Integer.class)
                .help("For positions suspected on homopolymer error if "
                + "total sum of frequencies with gap on this position "
                + "is [dominparamgen] times greater than without (or vice versa) "
                + "then error will be corrected. (Default: "
                + dominparamgen + ")");

        parser.addArgument("-dpp").dest(DOMINPOSTPROC_PARAMETER)
                .metavar("dominparampostpr")
                .setDefault(dominparampostpr)
                .type(Integer.class)
                .help("The same as [dominparamgen] but for pairwise "
                + "postprocessing of haplotypes with the help of "
                + "heuristics neigbor leaves of the phylogenetic tree."
                + "(Default: " + dominparampostpr + ")");

        try {
            Namespace n = parser.parseArgs(args);
            k = n.getInt(K_PARAMETER);
            nIter = n.getInt(I_PARAMETER);
            toFindHapl = n.getBoolean(ALIGN_PARAMETER);
            errorsseglen = n.getInt(L_PARAMETER);
            fl = (File) n.get(READS_PARAMETER);
            dominparamgen = n.getInt(DOMINGEN_PARAMETER);
            dominparampostpr = n.getInt(DOMINPOSTPROC_PARAMETER);
        } catch (ArgumentParserException e) {
            parser.handleError(e);
            System.exit(1);
        }

        int lt = 50;
        int mErPerc = 40;
        int toClust = 1;
        int dominparamonenucl = 5;
        int nucldiffparam = 1;
        int minNReads = 900;
        boolean toPrintStat = true;
        int maxz = 3;
        int kmin = k - 15;
        boolean toCalcHapl = true;
        boolean toDelUncor = false;
        boolean toPostprocessHeur = false; //???

        String dset_file = fl.getName();
        String outdir = "results" + "(" + dset_file + ")";
        DataSet ds = new DataSet(fl.getAbsolutePath());
        ds.setK(k);
        ds.setLenThr(lt);
        ds.setMaxAllErrorsPerc(mErPerc);
        ds.setFindErrorsSeglen(errorsseglen);

        Corrector cr = new Corrector(ds);
        cr.setMaxz(maxz);
        cr.setNIter(nIter);
        cr.setToClust(toClust);
        cr.setKmin(kmin);
        cr.setToRemoveAllUncorrect(toDelUncor);
        cr.setToFindHapl(toCalcHapl);
        cr.setToPrintStat(toPrintStat);
        cr.setToPostprocessHeur(toPostprocessHeur);
        cr.setMinNReads(minNReads);
        cr.run();
        ds = cr.CorrectedReads();

        String dset_file_name = fl.getAbsolutePath();
        ds.PrintCorrectedReads(fl.getAbsolutePath() + "_corrected.fas");

        // SECOND RUN
        maxz = 0;
        nIter = 0;
        toPrintStat = false;
        toDelUncor = true;
        toPostprocessHeur = false; //???
        dset_file_name = dset_file_name + "_corrected.fas";
        dset_file = dset_file_name;
        
        outdir = "results" + "(" + dset_file + ")";
        ds = new DataSet(dset_file);
        ds.setK(k);
        ds.setLenThr(lt);
        ds.setMaxAllErrorsPerc(mErPerc);
        ds.setFindErrorsSeglen(errorsseglen);
        minNReads = 1;


        cr = new Corrector(ds);
        cr.setMaxz(maxz);
        cr.setNIter(nIter);
        cr.setToClust(toClust);
        cr.setKmin(kmin);
        cr.setToRemoveAllUncorrect(toDelUncor);
        cr.setToFindHapl(toCalcHapl);
        cr.setToPrintStat(toPrintStat);
        cr.setToPostprocessHeur(toPostprocessHeur);
        cr.setMinNReads(minNReads);
        cr.run();
        ds = cr.CorrectedReads();

        ds.PrintCorrectedReads(dset_file_name + "_corrected.fas");
        ds.PrintHaplotypes(dset_file_name + "_haplotypes.fas");
        
        if (toFindHapl) {
            cr.postprocessHaplotypes(dset_file_name + "_haplotypes.fas", 15, 6.6, dominparampostpr);
            cr.printRevComp(dset_file_name + "_haplotypes.fas_postprocessed.fas");
            cr.postprocessHaplotypesPairwise(dset_file_name + "_haplotypes.fas_postprocessed.fas_RevComp.fas", 15, 6.6, dominparamonenucl, dominparamgen, nucldiffparam);
            cr.postprocessHaplotypes(dset_file_name + "_haplotypes.fas_postprocessed.fas_RevComp.fas_PostprocPair.fas", 15, 6.6, dominparampostpr);
            cr.postprocessHaplotypesPairwise(dset_file_name + "_haplotypes.fas_postprocessed.fas_RevComp.fas_PostprocPair.fas_postprocessed.fas", 15, 6.6, dominparamonenucl, dominparamgen, nucldiffparam);
        }
        System.out.println("Finished!");
        System.exit(0);
    }
}
