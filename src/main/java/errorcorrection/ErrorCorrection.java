package errorcorrection;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.security.CodeSource;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.argparse4j.*;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.*;

public class ErrorCorrection {

    private static final String K_PARAMETER = "K";
    private static final String I_PARAMETER = "I";
    private static final String L_PARAMETER = "L";
    private static final String READS_PARAMETER = "reads";
    private static final String DOMINGEN_PARAMETER = "domgen";
    private static final String DOMINPOSTPROC_PARAMETER = "dompostproc";
    private static final String CLUSTAL_PARAMETER = "clustalw";
    private static final String MUSCLE_PARAMETER = "muscle";
    static final String CLUSATL = "Clustal";
    static final String MUSCLE = "Muscle";
    static File env_path;
    static {
        try {
            CodeSource codeSource = ErrorCorrection.class.getProtectionDomain().getCodeSource();
            String path = codeSource.getLocation().toURI().getPath();
            env_path = new File(path).getParentFile();
        } catch (URISyntaxException ex) {
            Logger.getLogger(ErrorCorrection.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException, InterruptedException {
        ArgumentParser parser = ArgumentParsers.newArgumentParser("KEC.jar")
                .description("k-mer Error Correction algorithm input arguments parser.");

        int dominparamgen = 30;
        int dominparampostpr = 30; //30

        int k = 25;
        int nIter = 3;
        int errorsseglen = 0;
        File fl = null;
        String alignmethod = ""; // "Muscle" or "Clustal"

        MutuallyExclusiveGroup group = parser.addMutuallyExclusiveGroup("Align");

        group.addArgument("-cl","--clustalw").dest(CLUSTAL_PARAMETER)
                .action(Arguments.storeConst())
                .metavar("")
                .setConst(CLUSATL)
                .setDefault("")
                .type(String.class)
                .help("Enable using of CLustalW for multiple sequence "
                + "alignment for additional correction procedure (Default: "
                + "do not use)");
        group.addArgument("-mus","--muscle").dest(MUSCLE_PARAMETER)
                .action(Arguments.storeConst())
                .metavar("")
                .setConst(MUSCLE)
                .setDefault("")
                .type(String.class)
                .help("Enable using of Muscle for multiple sequence "
                + "alignment for additional correction procedure (Default: "
                + "do not use)");
        
        parser.addArgument(READS_PARAMETER)
                .metavar("ReadsFile")
                .help("File containing reads to be corrected "
                + " file with extension .fasta (.fas) "
                + "reads in fasta format")
                .type(File.class);

        parser.addArgument("-k").dest(K_PARAMETER)
                .metavar("SizeOfk-mer")
                .setDefault(k)
                .type(Integer.class)
                .help("Parameter k - the size of k-mers "
                + "Default: " + k + ")");

        parser.addArgument("-i").dest(I_PARAMETER)
                .metavar("Iterations")
                .setDefault(nIter)
                .type(Integer.class)
                .help("Number of iterations of the algorithm "
                + "(Default: " + nIter + ")");

        parser.addArgument("-l").dest(L_PARAMETER)
                .metavar("NumOfZeros")
                .setDefault(errorsseglen)
                .type(Integer.class)
                .help("Number of consecutive zeros for finding error threshold. "
                + "Using Poisson distribution if 0 (Default: " + errorsseglen + ")");

        parser.addArgument("-dg").dest(DOMINGEN_PARAMETER)
                .metavar("dominparamgen")
                .setDefault(dominparamgen)
                .type(Integer.class)
                .help("Parameter for haplotypes postprocessing using multiple alignment "
                + " For positions suspected on homopolymer error if "
                + "total sum of frequencies with gap on this position "
                + "is [dominparamgen] times greater than without (or vice versa) "
                + "then error will be corrected. (Default: "
                + dominparamgen + ")");

        parser.addArgument("-dpp").dest(DOMINPOSTPROC_PARAMETER)
                .metavar("dominparampostpr")
                .setDefault(dominparampostpr)
                .type(Integer.class)
                .help("The same as [dominparamgen] but for pairwise "
                + "postprocessing of haplotypes using "
                + "alignment of neigbor leaves of neighbor joining tree."
                + "(Default: " + dominparampostpr + ")");

        try {
            Namespace n = parser.parseArgs(args);
            k = n.getInt(K_PARAMETER);
            nIter = n.getInt(I_PARAMETER);
            String muscle = n.getString(MUSCLE_PARAMETER);
            String clustalw = n.getString(CLUSTAL_PARAMETER);
            alignmethod = "".equals(muscle) ? clustalw : muscle;
            errorsseglen = n.getInt(L_PARAMETER);
            fl = (File) n.get(READS_PARAMETER);
            dominparamgen = n.getInt(DOMINGEN_PARAMETER);
            dominparampostpr = n.getInt(DOMINPOSTPROC_PARAMETER);
        } catch (ArgumentParserException e) {
            parser.handleError(e);
            System.exit(1);
        }
        {
            System.out.println(alignmethod);
            int f = 100;
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
//		ds.setFreqThr(f);
            ds.setLenThr(lt);
            ds.setMaxAllErrorsPerc(mErPerc);
            ds.setFindErrorsSeglen(errorsseglen);

            Corrector cr = new Corrector(ds);
//                cr.printRevComp("Prep1-9.fas_corrected.fas_haplotypes.fas");
//                cr.postprocessHaplotypesPairwise("HOC_34_1b.fas_corrected.fas_haplotypes.fas_postprocessed.fas_RevComp.fas_PostprocPair.fas_postprocessed.fas_PostprocPair.fas", 15, 6.6, 3, 3, nucldiffparam);
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

            /*                cr = new Corrector(ds);
             cr.setMaxz(maxz);
             cr.setNIter(1);
             cr.setToClust(toClust);
             cr.setKmin(kmin);
             cr.setToRemoveAllUncorrect(false);
             cr.setToFindHapl(toCalcHapl);
             cr.setToPrintStat(toPrintStat);
             cr.run();
             ds = cr.CorrectedReads();*/


            String dset_file_name = fl.getAbsolutePath();
            ds.PrintCorrectedReads(fl.getAbsolutePath() + "_corrected.fas");
            /*                if (toFindHapl == 1)
             {
             ds.PrintHaplotypes(dset_file_name+"_haplotypes.fas");
             cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas",15,6.6, dominparampostpr);
             //                    cr.postprocessHaplotypes(dset_file_name+"_haplotypes.fas_postprocessed.fas",1,1);
             cr.printRevComp(dset_file_name +"_haplotypes.fas_postprocessed.fas");
             }
				
             System.out.println("Finished!");*/


            // SECOND RUN
            maxz = 0;
            nIter = 0;
            toPrintStat = false;
            toDelUncor = true;
            toPostprocessHeur = false; //???
            dset_file_name = dset_file_name + "_corrected.fas";
            dset_file = dset_file_name;
//                toFindHapl = 1;

            outdir = "results" + "(" + dset_file + ")";
            //System.out.println(dset_file);
            ds = new DataSet(dset_file);
            ds.setK(k);
//		ds.setFreqThr(f);
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
            if (!"".equals(alignmethod)) {
                /*                    dominparamgen = 5;
                 dominparampostpr = 5;
                 dominparamonenucl = 5;
                 nucldiffparam = 3;*/
                cr.postprocessHaplotypes(dset_file_name + "_haplotypes.fas", 15, 6.6, dominparampostpr,alignmethod);
//                cr.printRevComp(dset_file_name + "_haplotypes.fas_postprocessed.fas");
                cr.postprocessHaplotypesPairwise(dset_file_name + "_haplotypes.fas_postprocessed.fas", 15, 6.6, dominparamonenucl, dominparamgen, nucldiffparam,alignmethod);
                cr.postprocessHaplotypes(dset_file_name + "_haplotypes.fas_postprocessed.fas_PostprocPair.fas", 15, 6.6, dominparampostpr,alignmethod);
                cr.postprocessHaplotypesPairwise(dset_file_name + "_haplotypes.fas_postprocessed.fas_PostprocPair.fas_postprocessed.fas", 15, 6.6, dominparamonenucl, dominparamgen, nucldiffparam,alignmethod);
            }
//		cr.printRevComp(dset_file_name +"_haplotypes.fas_postprocessed.fas");

            System.out.println("Finished!");

        }


    }
}
