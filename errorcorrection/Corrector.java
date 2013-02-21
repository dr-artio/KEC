package errorcorrection;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Stack;
import java.util.StringTokenizer;

public class Corrector {
	DataSet ds;
	int kmin;
	int nIter;
	int maxz;
        int minNReads;
	boolean toClust;
	boolean toFindHapl;
        boolean toRemoveAllUncorrect;
        boolean toPrintStat;
        boolean toPostprocessHeur;
 	Corrector(DataSet d)
	{
		ds = d;
		kmin = ds.k - 10;
		nIter = 4;
		toClust = false;
		maxz = 1;
		toFindHapl = true;
                minNReads = 1000;
	}
        public Corrector()
        {

        }
	void setNIter(int j)
	{
		nIter = j;
	}
	void setMaxz(int i)
	{
		maxz = i;
	}
	void setKmin(int i)
	{
		kmin = i;
	}
	void setToFindHapl(boolean b)
	{
		toFindHapl = b;
	}
	void setToClust(int i)
	{
		if (i == 0)
			toClust = false;
		else
			toClust = true;
	}
        void setToRemoveAllUncorrect(boolean b)
        {
            toRemoveAllUncorrect = b;
        }
        void setMinNReads(int i)
        {
             minNReads = i;
        }
        void setToPrintStat(boolean b)
        {
            toPrintStat = b;
        }
        void setToPostprocessHeur(boolean b)
        {
            toPostprocessHeur = b;
        }
	void run() throws IOException
	{
		ArrayList<Integer> nErrors = new ArrayList<Integer>();
		ArrayList<Integer> maxReadsLen = new ArrayList<Integer>();
		ArrayList<Integer> minReadsLen = new ArrayList<Integer>();
		ArrayList<Double> AvReadsLen = new ArrayList<Double>();
		ArrayList<Integer> NReads = new ArrayList<Integer>();
		ArrayList<Integer> NUnReads = new ArrayList<Integer>();
		ArrayList<Integer> NUnBadReads = new ArrayList<Integer>();
                ArrayList<Integer> thresholds = new ArrayList<Integer>();
                ArrayList<Integer> readssizes = new ArrayList<Integer>();
                ArrayList<Integer> kmerssizes = new ArrayList<Integer>();
                ArrayList<Integer> allkmerssizes = new ArrayList<Integer>();
                ArrayList<Integer> usedmemory = new ArrayList<Integer>();
                ArrayList<Integer> Ncutted = new ArrayList<Integer>();
                Runtime rt = Runtime.getRuntime ();
                int mb = 1024*1024;

		
		ds.deleteStrangeNucl();
		calculations();
                if (toPrintStat)
                {
                    ds.PrintReadsStat("Statistics_" + ds.file_name);
//                    ds.printErrorsStat("Statistics_" + ds.file_name);
                }
		nErrors.add(ds.getNerrors());
		maxReadsLen.add(ds.maxReadLen());
		minReadsLen.add(ds.minReadLen());
		NReads.add(ds.getNreads());
		NUnReads.add(ds.getNuniquereads());
		NUnBadReads.add(ds.getNUniqueErrorReads());
		AvReadsLen.add(ds.AvReadLen());
                thresholds.add(ds.freqThr);
                readssizes.add(ds.reads.size());
                allkmerssizes.add(ds.allKmers.size());
                int kmerssz = 0;
                for (Read r : ds.reads)
                    kmerssz += r.kmers.size();
                kmerssizes.add(kmerssz);
                int ito = (int) (rt.totalMemory()/mb);
                int ifr = (int) (rt.freeMemory()/mb);
                usedmemory.add(ito - ifr);
                Ncutted.add(0);
		int i;
		
		if (!toClust)
		{
			for (int j = 0; j < nIter; j++)
			{
				i = ds.k;
				while (i >= kmin)
				{
					correctErrors(i);
					ds = new DataSet(ds);
					calculations();
					nErrors.add(ds.getNerrors());
					maxReadsLen.add(ds.maxReadLen());
					minReadsLen.add(ds.minReadLen());
					AvReadsLen.add(ds.AvReadLen());
					NReads.add(ds.getNreads());
					NUnReads.add(ds.getNuniquereads());
					NUnBadReads.add(ds.getNUniqueErrorReads());
                                        thresholds.add(ds.freqThr);
					i--;
				}
				cutTails();
				delLargeErrors();
				ds = new DataSet(ds);
				calculations();
				nErrors.add(ds.getNerrors());
				maxReadsLen.add(ds.maxReadLen());
				minReadsLen.add(ds.minReadLen());
				AvReadsLen.add(ds.AvReadLen());
				NReads.add(ds.getNreads());
				NUnReads.add(ds.getNuniquereads());
				NUnBadReads.add(ds.getNUniqueErrorReads());
                                thresholds.add(ds.freqThr);
			}
			
			correctErrors_heuristics();
			ds = new DataSet(ds);
			calculations();
			nErrors.add(ds.getNerrors());
			maxReadsLen.add(ds.maxReadLen());
			minReadsLen.add(ds.minReadLen());
			AvReadsLen.add(ds.AvReadLen());
			NReads.add(ds.getNreads());
			NUnReads.add(ds.getNuniquereads());
			NUnBadReads.add(ds.getNUniqueErrorReads());
                        thresholds.add(ds.freqThr);
                        thresholds.add(ds.freqThr);
                        readssizes.add(ds.reads.size());
                        allkmerssizes.add(ds.allKmers.size());
                        kmerssz = 0;
                        for (Read r : ds.reads)
                            kmerssz += r.kmers.size();
                        kmerssizes.add(kmerssz);
			
			
			for (int z = 0; z < maxz; z++)
			{
				for (int j = 0; j < nIter; j++)
				{
					i = ds.k;
					while (i >= kmin)
					{
						correctErrors(i);
						ds = new DataSet(ds);
						calculations();
						nErrors.add(ds.getNerrors());
						maxReadsLen.add(ds.maxReadLen());
						minReadsLen.add(ds.minReadLen());
						AvReadsLen.add(ds.AvReadLen());
						NReads.add(ds.getNreads());
						NUnReads.add(ds.getNuniquereads());
						NUnBadReads.add(ds.getNUniqueErrorReads());
                                                thresholds.add(ds.freqThr);
						i--;
					}
					
					correctErrors_heuristics();
					ds = new DataSet(ds);
					calculations();
					nErrors.add(ds.getNerrors());
					maxReadsLen.add(ds.maxReadLen());
					minReadsLen.add(ds.minReadLen());
					AvReadsLen.add(ds.AvReadLen());
					NReads.add(ds.getNreads());
					NUnReads.add(ds.getNuniquereads());
					NUnBadReads.add(ds.getNUniqueErrorReads());
                                        thresholds.add(ds.freqThr);
                                        readssizes.add(ds.reads.size());
                                        allkmerssizes.add(ds.allKmers.size());
                                        kmerssz = 0;
                                        for (Read r : ds.reads)
                                            kmerssz += r.kmers.size();
                                        kmerssizes.add(kmerssz);
					
					cutTails();
					delLargeErrors();
					ds = new DataSet(ds);
					calculations();
					nErrors.add(ds.getNerrors());
					maxReadsLen.add(ds.maxReadLen());
					minReadsLen.add(ds.minReadLen());
					AvReadsLen.add(ds.AvReadLen());
					NReads.add(ds.getNreads());
					NUnReads.add(ds.getNuniquereads());
					NUnBadReads.add(ds.getNUniqueErrorReads());
                                        thresholds.add(ds.freqThr);
                                        readssizes.add(ds.reads.size());
                                        allkmerssizes.add(ds.allKmers.size());
                                        kmerssz = 0;
                                        for (Read r : ds.reads)
                                            kmerssz += r.kmers.size();
                                        kmerssizes.add(kmerssz);
				}
				
			}
		}
		else
		{
//                        String refr = "GCATGGGATATGATGATGAACTGGAGTCCAACGACCACCTTACTCCTCGCCCAGGTTATGAGGATCCCAGGTACTCTGGTAGATTTACTCGCTGGAGGCCACTGGGGTGTCCTCGTGGGAGTGGCCTATTTCAGTATGCAAGCCAACTGGGCCAAAGTCATCTTGGTCCTATTCCTTTTTGCAGGGGTTGACGCCAAGACCACCACAACTGGGTCTGCGGCGGCCAGGGGAGTCAGCCGCGTCACTGGGTTCTTTGCCCCCGGGCCTAGCCAAAACTTGCAGCTCATTAACACCAACGGGAGCTGGCA";
//                        int d = findClosest(refr,15,6.6);
//                        System.out.println();
                        
			for (int j = 0; j < nIter; j++)
			{
                            
                                
                                i = ds.k;
				while (i >= kmin)
				{
					correctErrors(i);
					i--;
				}
				int ncut = cutTails();
				delLargeErrors();
				ds = new DataSet(ds);
//                                d = findClosest(refr,15,6.6);
//                                int ncut = 0;
                                System.gc();
				calculations();
				nErrors.add(ds.getNerrors());
				maxReadsLen.add(ds.maxReadLen());
				minReadsLen.add(ds.minReadLen());
				AvReadsLen.add(ds.AvReadLen());
				NReads.add(ds.getNreads());
				NUnReads.add(ds.getNuniquereads());
				NUnBadReads.add(ds.getNUniqueErrorReads());
                                thresholds.add(ds.freqThr);
                                readssizes.add(ds.reads.size());
                                allkmerssizes.add(ds.allKmers.size());
                                kmerssz = 0;
                                for (Read r : ds.reads)
                                    kmerssz += r.kmers.size();
                                kmerssizes.add(kmerssz);
                                ito = (int) (rt.totalMemory()/mb);
                                ifr = (int) (rt.freeMemory()/mb);
                                usedmemory.add(ito - ifr);
                                Ncutted.add(ncut);
                                if (ds.getNreads() < minNReads)
                                    break;

			}

                        if (toPostprocessHeur)
                        {
                            correctErrors_heuristics();
                            ds = new DataSet(ds);
                            System.gc();
                            calculations();
                            nErrors.add(ds.getNerrors());
                            maxReadsLen.add(ds.maxReadLen());
                            minReadsLen.add(ds.minReadLen());
                            AvReadsLen.add(ds.AvReadLen());
                            NReads.add(ds.getNreads());
                            NUnReads.add(ds.getNuniquereads());
                            NUnBadReads.add(ds.getNUniqueErrorReads());
                            thresholds.add(ds.freqThr);
                            readssizes.add(ds.reads.size());
                            allkmerssizes.add(ds.allKmers.size());
                            kmerssz = 0;
                            for (Read r : ds.reads)
                                kmerssz += r.kmers.size();
                            kmerssizes.add(kmerssz);
                            ito = (int) (rt.totalMemory()/mb);
                            ifr = (int) (rt.freeMemory()/mb);
                            usedmemory.add(ito - ifr);
			
			
                            for (int z = 0; z < maxz; z++)
                            {
				for (int j = 0; j < nIter; j++)
				{
					i = ds.k;
					while (i >= kmin)
					{
						correctErrors(i);
						i--;
					}				
					correctErrors_heuristics();
					cutTails();
					delLargeErrors();
					ds = new DataSet(ds);
                                        System.gc();
					calculations();
					nErrors.add(ds.getNerrors());
					maxReadsLen.add(ds.maxReadLen());
					minReadsLen.add(ds.minReadLen());
					AvReadsLen.add(ds.AvReadLen());
					NReads.add(ds.getNreads());
					NUnReads.add(ds.getNuniquereads());
					NUnBadReads.add(ds.getNUniqueErrorReads());
                                        thresholds.add(ds.freqThr);
                                        readssizes.add(ds.reads.size());
                                        allkmerssizes.add(ds.allKmers.size());
                                        kmerssz = 0;
                                        for (Read r : ds.reads)
                                            kmerssz += r.kmers.size();
                                        kmerssizes.add(kmerssz);
                                        ito = (int) (rt.totalMemory()/mb);
                                        ifr = (int) (rt.freeMemory()/mb);
                                        usedmemory.add(ito - ifr);
				}				
                            }
                            for (int j = 0; j < nIter; j++)
                            {
                            	i = ds.k;
				while (i >= kmin)
				{
					correctErrors(i);
					i--;
				}
				cutTails();
				delLargeErrors();
				ds = new DataSet(ds);
                                System.gc();
				calculations();
				nErrors.add(ds.getNerrors());
				maxReadsLen.add(ds.maxReadLen());
				minReadsLen.add(ds.minReadLen());
				AvReadsLen.add(ds.AvReadLen());
				NReads.add(ds.getNreads());
				NUnReads.add(ds.getNuniquereads());
				NUnBadReads.add(ds.getNUniqueErrorReads());
                                thresholds.add(ds.freqThr);
                                readssizes.add(ds.reads.size());
                                allkmerssizes.add(ds.allKmers.size());
                                kmerssz = 0;
                                for (Read r : ds.reads)
                                    kmerssz += r.kmers.size();
                                kmerssizes.add(kmerssz);
                                ito = (int) (rt.totalMemory()/mb);
                                ifr = (int) (rt.freeMemory()/mb);
                                usedmemory.add(ito - ifr);
                            }
                    }
		}
		
		
//		ds.PrintReadsStat("Statistics_after_" + ds.file_name);
//		ds.printErrorsStat("Statistics_after_" + ds.file_name);


                if ((toRemoveAllUncorrect)&&(ds.getNreads() > minNReads))
                    deleteWholeUncorrectible();

                ds = new DataSet(ds);
                
//                String refr = "GCATGGGATATGATGATGAACTGGAGTCCAACGACCACCTTACTCCTCGCCCAGGTTATGAGGATCCCAGGTACTCTGGTAGATTTACTCGCTGGAGGCCACTGGGGTGTCCTCGTGGGAGTGGCCTATTTCAGTATGCAAGCCAACTGGGCCAAAGTCATCTTGGTCCTATTCCTTTTTGCAGGGGTTGACGCCAAGACCACCACAACTGGGTCTGCGGCGGCCAGGGGAGTCAGCCGCGTCACTGGGTTCTTTGCCCCCGGGCCTAGCCAAAACTTGCAGCTCATTAACACCAACGGGAGCTGGCA";
//                int d = findClosest(refr,15,6.6);
                
                System.gc();
		calculations();
		nErrors.add(ds.getNerrors());
		maxReadsLen.add(ds.maxReadLen());
		minReadsLen.add(ds.minReadLen());
		AvReadsLen.add(ds.AvReadLen());
		NReads.add(ds.getNreads());
		NUnReads.add(ds.getNuniquereads());
		NUnBadReads.add(ds.getNUniqueErrorReads());
                thresholds.add(ds.freqThr);
                readssizes.add(ds.reads.size());
                allkmerssizes.add(ds.allKmers.size());
                kmerssz = 0;
                for (Read r : ds.reads)
                    kmerssz += r.kmers.size();
                kmerssizes.add(kmerssz);
                ito = (int) (rt.totalMemory()/mb);
                ifr = (int) (rt.freeMemory()/mb);
                usedmemory.add(ito - ifr);

                System.out.println("Thresholds: " + thresholds);
                System.out.println("Errors: " + nErrors);
		System.out.println("MaxReadslen: " + maxReadsLen);
		System.out.println("MinReadslen: " + minReadsLen);
		System.out.println("AvReadslen: " + AvReadsLen);
		System.out.println("NReads: " + NReads);
		System.out.println("NUnReads: " + NUnReads);
		System.out.println("NUnBadReads: " + NUnBadReads);
                System.out.println("NCutted: " + Ncutted);

		if (toFindHapl)
			ds.findHaplotypes();
		System.out.println("Haplotypes: " + ds.haplotypes.size());
                
 //               d = findClosestHapl(refr,15,6.6);

                FileWriter fw = new FileWriter(ds.file_name + "_log.txt");
                fw.write("Thresholds: " + thresholds + "\n");
                fw.write("Errors: " + nErrors + "\n");
                fw.write("MaxReadslen: " + maxReadsLen+ "\n");
		fw.write("MinReadslen: " + minReadsLen+ "\n");
		fw.write("AvReadslen: " + AvReadsLen+ "\n");
		fw.write("NReads: " + NReads+ "\n");
		fw.write("NUnReads: " + NUnReads+ "\n");
		fw.write("NUnBadReads: " + NUnBadReads+ "\n");
                fw.write("Haplotypes: " + ds.haplotypes.size()+ "\n");
                fw.write("ReadsSizes: " + readssizes+ "\n");
                fw.write("AllKmersSizes: " + allkmerssizes+ "\n");
                fw.write("KmersSizes: " + kmerssizes+ "\n");
                fw.write("UsedMemory: " + usedmemory+ "\n");
                fw.write("NCutted: " + Ncutted+ "\n");
                fw.close();
	}
	DataSet CorrectedReads()
	{
		return ds;
	}
	void correctErrors(int length)
	{
                int count = 0;
		for (Read r : ds.badReads)
                {
                        count++;
			for (ErrorRegion er : r.errorRegions)
				if (er.length == length)
				{
					if ((er.length == ds.k)&&(er.begin!= 0)&&(er.end + ds.k != r.getLength()))
					{
						String s = r.nucl.substring(er.begin, er.end);
						char c = r.nucl.charAt(er.end);
						char repl = c;
						int maxKcount = 0;
						for (Kmer_general km : ds.allKmers)
							if ((km.kCount > ds.freqThr)&&(km.nucl.startsWith(s)))
								if (km.nucl.charAt(ds.k-1)!= c)
									if (km.kCount > maxKcount)
									{
										maxKcount = km.kCount;
										repl = km.nucl.charAt(ds.k-1);
									}
/*						if (maxKcount > ds.freqThr)
							r.corrections[er.end] = new Correction('R',""+repl);
						else
							r.corrections[er.end] = new Correction('D'); */
                                                if (maxKcount > ds.freqThr)
                                                    if (repl != r.getNucl(er.end + 1))
							r.corrections[er.end] = new Correction('R',""+repl);
                                                    else
                                                    {
                                                        s = r.getNucl(er.begin + 1, er.end -1) + repl + repl;
                                                        char q = 'D';
                                                        for (Kmer_general km : ds.allKmers)
                                                            if (km.nucl.equalsIgnoreCase(s)&& (km.kCount > ds.freqThr))
                                                            {
                                                                q = 'R';
                                                                break;
                                                            }
                                                        if (q == 'D')
                                                            r.corrections[er.end] = new Correction('D');
                                                        if (q == 'R')
                                                            r.corrections[er.end] = new Correction('R',""+repl);
                                                    }
						else
							r.corrections[er.end] = new Correction('D');
					}
					if ((er.length == ds.k-1)&&(er.begin!= 0)&&(er.end + ds.k < r.getLength()-1))
					{
                                                    String s = r.getNucl(er.begin, er.end);
                                                    char let[] = new char[4];
                                                    let[0] = 'A';
                                                    let[1] = 'T';
                                                    let[2] = 'G';
                                                    let[3] = 'C';
                                                    int bestKcount = 0;
                                                    char bestins = 'W';
                                                    for (Kmer_general km : ds.allKmers)
							for (int i = 0; i < 4; i++)
								if ((km.kCount > ds.freqThr)&&(km.nucl.contentEquals(s+let[i])))
								{
									if (km.kCount > bestKcount)
									{
										bestKcount = km.kCount;
										bestins = let[i];
									}
								}
                                                
						if (r.getNucl(er.end) != r.getNucl(er.end+1))  
                                                {
                                            
							if (bestins != 'W')
								r.corrections[er.end] = new Correction('I',bestins);
                                                }
						else
                                                {
                                                        if ((bestins != r.getNucl(er.end + 2)) && (bestins != 'W'))
                                                            r.corrections[er.end] = new Correction('I',bestins);
                                                        else
                                                            r.corrections[er.end] = new Correction('D');
                                                }
					}
					if (er.length <= ds.k - 2)
					{
						int hom_ins = ds.k + 1 - er.length;
						int hom_del = ds.k - 1 - er.length;
						int l = 0;
						int q = 0;
						int i = er.end+1;
						int j = er.end + 2;
						while((r.getNucl(i) == r.getNucl(er.end))&&(i < r.getLength()))
						{
							l++;
							i++;
						}
						if (r.getNucl(er.end)!= r.getNucl(er.end + 1))
							while((j < r.getLength()) && (r.getNucl(j) == r.getNucl(er.end+1)))
							{
								q++;
								j++;
							}
						int bestKcount = 0;
						if (q == hom_del - 1)
						{
							String s = r.getNucl(er.begin, er.end+1+q);
							s = s + r.getNucl(er.end+1);
							
							for (Kmer_general km : ds.allKmers)
								if ((km.kCount > ds.freqThr)&&(km.nucl.contentEquals(s)))
								{
									if (km.kCount > bestKcount)
										bestKcount = km.kCount;
								}
						}
						if (l == hom_ins -1)
							r.corrections[er.end] = new Correction('D');
						if (q==hom_del - 1)
							if (bestKcount > 0)
								r.corrections[er.end+1+q] = new Correction('I',r.getNucl(er.end+1));
					}	
				}
            }
            DynamicOut.finishSteps();
	}
	void correctErrors_heuristics()
	{
		for (Read r : ds.badReads)
			for (ErrorRegion er : r.errorRegions)
				if (er.length > ds.k)
				{
					boolean corr = false;
					if (r.getNucl(er.begin + ds.k - 2) == r.getNucl(er.begin + ds.k - 1))
					{
						r.corrections[er.begin + ds.k - 1] = new Correction('D');
						continue;
					}
					if (r.getNucl(er.end) == r.getNucl(er.end+1))
					{
						r.corrections[er.end] = new Correction('D');
						continue;
					}
					if (corr)
						continue;
					if (er.end + 2 < r.getLength())
						if (r.getNucl(er.end) != r.getNucl(er.end+1))
						{							
							int j = er.end + 2;
							while((r.getNucl(j) == r.getNucl(er.end+1))&&(j < r.getLength()))
								j++;
							int bestKcount = 0;
							String s = r.getNucl(j-ds.k+1, j-1);
							s = s + r.getNucl(er.end+1);
								
							for (Kmer_general km : ds.allKmers)
								if ((km.kCount > ds.freqThr)&&(km.nucl.contentEquals(s)))
								{
									if (km.kCount > bestKcount)
										bestKcount = km.kCount;
								}
							if (bestKcount > 0)
                                                        {
								r.corrections[er.end+1] = new Correction('I',r.getNucl(er.end+1));
                                                                continue;
                                                        }
//							r.corrections[er.end] = new Correction('I', r.getNucl(er.end + 1));
						}
//					if ((r.getNucl(er.begin + ds.k - 2) != r.getNucl(er.begin + ds.k - 1))&&(r.getNucl(er.begin + ds.k - 1) == r.getNucl(er.begin + ds.k)))
//							r.corrections[er.begin + ds.k - 2] = new Correction('I', r.getNucl(er.begin + ds.k - 1));
					if (corr)
						continue;
					
					String s = r.nucl.substring(er.begin, er.begin + ds.k - 1);
					char c = r.nucl.charAt(er.begin + ds.k - 1);
					char repl = c;
					int maxKcount = 0;
					for (Kmer_general km : ds.allKmers)
						if ((km.kCount > ds.freqThr)&&(km.nucl.startsWith(s)))
							if (km.nucl.charAt(ds.k-1)!= c)
								if (km.kCount > maxKcount)
								{
									maxKcount = km.kCount;
									repl = km.nucl.charAt(ds.k-1);
								}
					if (maxKcount > ds.freqThr)
						r.corrections[er.begin + ds.k - 1] = new Correction('R',""+repl);
					else
						r.corrections[er.begin + ds.k - 1] = new Correction('D');		
					
					s = r.nucl.substring(er.end-ds.k+1, er.end);
					c = r.nucl.charAt(er.end);
					repl = c;
					maxKcount = 0;
					for (Kmer_general km : ds.allKmers)
						if ((km.kCount > ds.freqThr)&&(km.nucl.startsWith(s)))
							if (km.nucl.charAt(ds.k-1)!= c)
								if (km.kCount > maxKcount)
								{
									maxKcount = km.kCount;
									repl = km.nucl.charAt(ds.k-1);
								}
					if (maxKcount > ds.freqThr)
						r.corrections[er.end] = new Correction('R',""+repl);
					else
						r.corrections[er.end] = new Correction('D');	
				}
	}
	int cutTails()
	{
                if (ds.getNreads() < minNReads)
                    return 0;
		int count = 0;
		for (Read r : ds.badReads)
		{
	//		System.out.println("Deleting tails: " + r.name +"\n" );
			ArrayList<ErrorRegion> corrected = new ArrayList<ErrorRegion>();
			for (ErrorRegion er : r.errorRegions)
			{
				if (er.begin == 0)
				{
					for (int i = er.begin; i <= er.end; i++)
						r.corrections[i] = new Correction('D');
					count++;
				}
				
				if ((er.end + ds.k  == r.getLength())&&(er.begin != 0))
				{
					for (int i = er.begin+ds.k-1; i < r.getLength(); i++)
						r.corrections[i] = new Correction('D');
					count++;
				}
			}
		}
		System.out.println("Tails cutted: " + count);
                return count;
	}
	void delLargeErrors()
	{
                if (ds.getNreads() < minNReads)
                    return;
		int count = 0;
                for (Read r : ds.badReads)
		{
			double l = 0;
			int lr = r.getLength();
			for (ErrorRegion er : r.errorRegions)
				l+=er.length;
			double d =  100.0*(double)(l/lr);
			if (d > ds.maxAllErrorsPerc)
				count++;
		}
                if (ds.getNreads() - count < minNReads)
                    return;
                count = 0;
		for (Read r : ds.badReads)
		{
//                    if (r.name.equalsIgnoreCase("read385"))
//                        System.out.println();
			double l = 0;
			int lr = r.getLength();
			for (ErrorRegion er : r.errorRegions)
				l+=er.length;
			double d =  100.0*(double)(l/lr);
//			if (d > r.getLength()*(ds.maxAllErrorsPerc/100))
			if (d > ds.maxAllErrorsPerc)
			{
				ds.reads.remove(r);
				count++;
			}
		}
		System.out.println("Reads with large errors deleted: " + count);
	}
	void calculations() throws IOException
	{
		ds.calculateKMersAndKCounts();

                if (ds.freqThr == -1)
                {
                    if (ds.finderrorsseglen == 0)
                        ds.findFreqThresholdPoisson();
                    else
                        ds.findFreqThreshold();
                }

		if (toClust)
		{
			ds.clustering();
			ds.clear();		
		}
		ds.findErrors();
	}
        void removeUncorrectible()
        {
            ArrayList<Read> newReads = new ArrayList<Read>();
            for (Read r : ds.badReads)
            {
                for (ErrorRegion er : r.errorRegions)
                {
                    if (er.length <= ds.k)
                        r.replaceNucl(er.end, 'E');
                    if ((er.length > ds.k)&&(er.end < 2*ds.k))
                    {
                        r.replaceNucl(er.begin + ds.k - 1, 'E');
                        r.replaceNucl(er.end, 'E');
                    }
                    if (er.length >= 2*ds.k)
                        for (int i = er.begin+ds.k-1; i <= er.end;i++)
                            r.replaceNucl(i, 'E');
                }
                StringTokenizer st = new StringTokenizer(r.nucl,"E");
                while (st.hasMoreTokens())
                {
                    String s = st.nextToken();
                    if (s.length() > ds.lengthThr)
                        newReads.add(new Read(s,r.frequency));
                }
            }
            ds.reads.removeAll(ds.badReads);
            ds.reads.addAll(newReads);
        }
        void deleteWholeUncorrectible()
        {
            ds.reads.removeAll(ds.badReads);
        }
        void EDAR() throws IOException
        {
                
                ArrayList<Integer> nErrors = new ArrayList<Integer>();
		ArrayList<Integer> maxReadsLen = new ArrayList<Integer>();
		ArrayList<Integer> minReadsLen = new ArrayList<Integer>();
		ArrayList<Double> AvReadsLen = new ArrayList<Double>();
		ArrayList<Integer> NReads = new ArrayList<Integer>();
		ArrayList<Integer> NUnReads = new ArrayList<Integer>();
		ArrayList<Integer> NUnBadReads = new ArrayList<Integer>();

		ds.deleteStrangeNucl();
		calculations();
                removeUncorrectible();
                ds = new DataSet(ds);
                calculations();
		nErrors.add(ds.getNerrors());
		maxReadsLen.add(ds.maxReadLen());
		minReadsLen.add(ds.minReadLen());
		AvReadsLen.add(ds.AvReadLen());
		NReads.add(ds.getNreads());
		NUnReads.add(ds.getNuniquereads());
		NUnBadReads.add(ds.getNUniqueErrorReads());

		System.out.println("Errors: " + nErrors);
		System.out.println("MaxReadslen: " + maxReadsLen);
		System.out.println("MinReadslen: " + minReadsLen);
		System.out.println("AvReadslen: " + AvReadsLen);
		System.out.println("NReads: " + NReads);
		System.out.println("NUnReads: " + NUnReads);
		System.out.println("NUnBadReads: " + NUnBadReads);
		if (toFindHapl)
			ds.findHaplotypes();
		System.out.println("Haplotypes: " + ds.haplotypes.size());
        }
        void postprocessHaplotypes(String addr, double gapop, double gapext, int dominparampostpr) throws IOException, InterruptedException
        {
                File fl = new File(addr);
                if (!fl.exists())
                {
                    System.out.println("No such file!");
                    return;
                }
                if ((ds!=null)&&(ds.haplotypes.size() == 1))
                {
                    System.out.println("No postprocessing needed");
                    File f1 = new File(addr);
                    File f2 = new File(addr + "_postprocessed.fas");
                    f1.renameTo(f2);
                    return;
                }
                DataSet dats = new DataSet(addr);
                if (dats.getNreads() == 1)
                {
                    System.out.println("No postprocessing needed");
                    return;
                }

                System.out.println("Postprocessing haplotypes \n");
                Runtime run=Runtime.getRuntime();
                Process p=null;
                String s = "ClustalW2//clustalw2 -INFILE=" +addr + " -OUTFILE=" + addr+"_allign.fas";
                s+= " -OUTPUT=FASTA -DNAMATRIX=IUB -GAPOPEN=" + gapop +" -GAPEXT=" + gapext + " -TYPE=DNA -PWDNAMATRIX=IUB -PWGAPOPEN=" + gapop+ " -PWGAPEXT=" + gapext;

                p = run.exec(s);
                BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));

                System.out.println("Here is the standard output of the command:\n");
                while ((s = stdInput.readLine()) != null) {
                    System.out.println(s);
                }
                p.waitFor();
                if (p.exitValue()!= 0)
                {
                    System.out.println("Error in allgnment program");
                }
                stdInput.close();
                
                DataSet hap = new DataSet(addr+"_allign.fas",'c');

                int len = 0;
                for (Read r : hap.reads)
                {
                    StringTokenizer st = new StringTokenizer(r.name,"|");
                    st.nextToken();
                    double d = 10000*Double.parseDouble(st.nextToken());
                    r.frequency = (int) d;
                    len = r.getLength();
                }
                for (int i = 0; i < len; i++)
                {
                    HashMap<Character,Integer> hm = new HashMap<Character,Integer>();
                    int count = 0;
                    for (Read r : hap.reads)
                    {
                        char c = r.getNucl(i);
                        if (hm.containsKey(c))
                            hm.put(c, hm.get(c) + r.frequency);
                        else
                        {
                            hm.put(c, r.frequency);
                            count++;
                        }
                    }
                    if ((count == 1)&&(hm.containsKey('-')))
                        for (Read r : hap.reads)
                            r.corrections[i] = new Correction('D');
                    if (count != 2)
                        continue;
                    if (hm.containsKey('-'))
                    {
                        char nongap = 'z';
                        int nongap_fr = 0;
                        int gap_fr = 0;
                        Iterator ir = hm.entrySet().iterator();
                        while (ir.hasNext())
                        {
                            Map.Entry me = (Map.Entry) ir.next();
                            if ((Character) me.getKey() == '-')
                                gap_fr = (Integer) me.getValue();
                            else
                            {
                                nongap_fr = (Integer) me.getValue();
                                nongap = (Character) me.getKey();
                            }
                        }
                        boolean homopError = true;
                        for (Read r : hap.reads)
                        {
                            if (i == 0)
                            {
                                homopError = false;
                                break;
                            }
                            if (r.getNucl(i-1) != nongap)
                            {
                                if (r.getNucl(i-1) != '-')
                                {
                                    homopError = false;
                                    break;
                                }
                            }
                        }
                        if (homopError)
                        {
                            if (gap_fr > dominparampostpr * nongap_fr)
                                for (Read r : hap.reads)
                                    r.corrections[i] = new Correction('D');
                            if (gap_fr * dominparampostpr <= nongap_fr)
                                for (Read r : hap.reads)
                                    r.corrections[i] = new Correction('R',nongap);
                        }
                        else
                        {
                            if (i == len -1)
                                continue;
                             homopError = true;
                             for (Read r : hap.reads)
                                if (r.getNucl(i+1) != nongap)
                                {
                                    homopError = false;
                                    break;
                                }
                             if (homopError)
                             {
                                if (gap_fr > nongap_fr * dominparampostpr)
                                    for (Read r : hap.reads)
                                        r.corrections[i] = new Correction('D');
                                if (gap_fr * dominparampostpr <= nongap_fr)
                                    for (Read r : hap.reads)
                                        r.corrections[i] = new Correction('R',nongap);
                            }
                        }
                    }

                }
              DataSet corhap = new DataSet(hap);
              corhap.findHaplotypes();
              StringTokenizer st1 = new StringTokenizer(addr,"_");
              String tag = st1.nextToken();
              corhap.PrintHaplotypesWithNameTag(addr+"_postprocessed.fas", tag);
//              corhap.PrintHaplotypes(addr+"_postprocessed.fas");
        }
         void postprocessHaplotypesSHORAH(String addr) throws IOException, InterruptedException
        {
                System.out.println("Postprocessing haplotypes \n");
                DataSet hap = new DataSet(addr,'c');

                int len = 0;
                for (Read r : hap.reads)
                {
                    StringTokenizer st = new StringTokenizer(r.name,"_");
                    st.nextToken();
                    double d = 10000*Double.parseDouble(st.nextToken());
                    r.frequency = (int) d;
                    len = r.getLength();
                }
                for (int i = 1; i < len; i++)
                {
                    HashMap<Character,Integer> hm = new HashMap<Character,Integer>();
                    int count = 0;
                    for (Read r : hap.reads)
                    {
                        char c = r.getNucl(i);
                        if (hm.containsKey(c))
                            hm.put(c, hm.get(c) + r.frequency);
                        else
                        {
                            hm.put(c, r.frequency);
                            count++;
                        }
                    }
                    if ((count == 1)&&(hm.containsKey('-')))
                        for (Read r : hap.reads)
                            r.corrections[i] = new Correction('D');
                    if (count != 2)
                        continue;
                    if (hm.containsKey('-'))
                    {
                        char nongap = 'z';
                        int nongap_fr = 0;
                        int gap_fr = 0;
                        Iterator ir = hm.entrySet().iterator();
                        while (ir.hasNext())
                        {
                            Map.Entry me = (Map.Entry) ir.next();
                            if ((Character) me.getKey() == '-')
                                gap_fr = (Integer) me.getValue();
                            else
                            {
                                nongap_fr = (Integer) me.getValue();
                                nongap = (Character) me.getKey();
                            }
                        }
                        boolean homopError = true;
                        for (Read r : hap.reads)
                            if (r.getNucl(i-1) != nongap)
                            {
//                                if (r.getNucl(i-1) != '-')
                                {
                                    homopError = false;
                                    break;
                                }
                            }
                        if (homopError)
                        {
                            if (gap_fr > nongap_fr)
                                for (Read r : hap.reads)
                                    r.corrections[i] = new Correction('D');
                            else
                                for (Read r : hap.reads)
                                    r.corrections[i] = new Correction('R',nongap);
                        }
                        else
                        {
                            if (i == len -1)
                                continue;
                             homopError = true;
                             for (Read r : hap.reads)
                                if (r.getNucl(i+1) != nongap)
                                {
                                    homopError = false;
                                    break;
                                }
                             if (homopError)
                             {
                                if (gap_fr > nongap_fr)
                                    for (Read r : hap.reads)
                                        r.corrections[i] = new Correction('D');
                                else
                                    for (Read r : hap.reads)
                                        r.corrections[i] = new Correction('R',nongap);
                            }
                        }
                    }
                }
              DataSet corhap = new DataSet(hap);
              corhap.findHaplotypes();
              corhap.PrintHaplotypes(addr+"_postprocessed.fas");
        }
        void delCompGapsSHORAH(String addr) throws IOException
        {
                System.out.println("Postprocessing haplotypes \n");
                DataSet hap = new DataSet(addr,'c');

                int len = 0;
                for (Read r : hap.reads)
                {
                    StringTokenizer st = new StringTokenizer(r.name,"_");
                    st.nextToken();
                    double d = 10000*Double.parseDouble(st.nextToken());
                    r.frequency = (int) d;
                    len = r.getLength();
                }
                for (int i = 1; i < len; i++)
                {
                    HashMap<Character,Integer> hm = new HashMap<Character,Integer>();
                    int count = 0;
                    for (Read r : hap.reads)
                    {
                        char c = r.getNucl(i);
                        if (hm.containsKey(c))
                            hm.put(c, hm.get(c) + r.frequency);
                        else
                        {
                            hm.put(c, r.frequency);
                            count++;
                        }
                    }
                    if ((count == 1)&&(hm.containsKey('-')))
                        for (Read r : hap.reads)
                            r.corrections[i] = new Correction('D');
                }
              hap.PrintCorrectedReads(addr+"delgap.fas");
        }
        void printRevComp(String addr) throws IOException
        {
                File f = new File(addr);
                if (!f.exists())
                {
                    System.out.println("No file to reverse");
                    return;
                }
                FileReader fr = new FileReader(addr);
		BufferedReader br = new BufferedReader(fr);
                FileWriter fw = new FileWriter(addr+"_RevComp.fas");
		
		String s = "";
                String nucl = "";
		while (s!= null)
		{
			if (s.length() == 0)
			{
				s = br.readLine();
				continue;
			}
			if(s.charAt(0) == '>')
			{
                            for (int i = nucl.length()-1; i>=0; i--)
                            {
                                if (nucl.charAt(i) == 'A')
                                    fw.write("T");
                                if (nucl.charAt(i) == 'T')
                                    fw.write("A");
                                if (nucl.charAt(i) == 'G')
                                    fw.write("C");
                                if (nucl.charAt(i) == 'C')
                                    fw.write("G");
                            }
                            fw.write("\n");
                            nucl = "";
                            fw.write(s+"\n");

			}
			else
				nucl+=s.toUpperCase();
			s = br.readLine();
		}
                for (int i = nucl.length()-1; i>=0; i--)
                {
                    if (nucl.charAt(i) == 'A')
                           fw.write("T");
                    if (nucl.charAt(i) == 'T')
                           fw.write("A");
                    if (nucl.charAt(i) == 'G')
                           fw.write("C");
                    if (nucl.charAt(i) == 'C')
                           fw.write("G");
                }
                fw.close();
    }
        void FindMaxSeq(String addr) throws IOException
        {
            DataSet ds = new DataSet(addr);
            ds.findHaplotypes();
            ds.PrintHaplotypes(addr+"MaxSeq.fas");
        }
       public void postprocessHaplotypesPairwise(String addr, double gapop, double gapext, int dominparamonenucl, int dominparamgen, int nucldiffparam) throws IOException, InterruptedException
       {
           int dominparamgenins = 100;
            File fl = new File(addr);
                if (!fl.exists())
                {
                    System.out.println("No such file!");
                    return;
                }
            DataSet hap = new DataSet(addr,'c');
            if (hap.reads.size() == 1)
            {
                File f1 = new File(addr);
                File f2 = new File(addr + "_PostprocPair.fas");
                f1.renameTo(f2);
                return;
            }

            int niterpostp = 20;
            for (int iter = 0; iter < niterpostp; iter++)
            {
                System.out.println("Postprocessing haplotypes pairwise \n");
                Runtime run=Runtime.getRuntime();
                Process p=null;
                String s = "ClustalW2//clustalw2 -INFILE=" +addr + " -OUTFILE=" + addr+"_allign.fas";
                s+= " -OUTPUT=FASTA -DNAMATRIX=IUB -GAPOPEN=" + gapop +" -GAPEXT=" + gapext + " -TYPE=DNA -PWDNAMATRIX=IUB -PWGAPOPEN=" + gapop+ " -PWGAPEXT=" + gapext;

                p = run.exec(s);
                BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));

                System.out.println("Here is the standard output of the command:\n");
                while ((s = stdInput.readLine()) != null) {
                    System.out.println(s);
                }
                p.waitFor();
                if (p.exitValue()!= 0)
                {
                    System.out.println("Error in allgnment program");
                }
                stdInput.close();
                int i = addr.length() - 1;
                while(addr.charAt(i) != '.')
                    i--;
                FileReader ftree = new FileReader(addr.substring(0, i) + ".dnd");
                BufferedReader brtree = new BufferedReader(ftree);
                for (Read r : hap.reads)
                {
                    StringTokenizer st = new StringTokenizer(r.name,"|");
                    st.nextToken();
                    double d = 10000*Double.parseDouble(st.nextToken());
                    r.frequency = (int) d;
                }
                s = brtree.readLine();
                String treestr = "";
                while (s!= null)
                {
                    treestr += s;
                    s = brtree.readLine();
                }
                brtree.close();
                HashSet<ArrayList<Read>> closepairs = new HashSet<ArrayList<Read>>();
                Stack<Character> st = new Stack<Character>();
                for (i = 0; i < treestr.length(); i++)
                {
                    st.push(treestr.charAt(i));
                    if (treestr.charAt(i) == ')')
                    {
                        String betw = "";
                        char c = ' ';
                        while (c != '(')
                        {
                            c = st.pop();
                            betw = c + betw;
                        }
                        ArrayList<String> readsbetw = new ArrayList<String>();
                        StringTokenizer stbetw = new StringTokenizer(betw,"()_,");
                        while (stbetw.hasMoreTokens())
                        {
                            String s1 = stbetw.nextToken();
                            if (s1.startsWith("read"))
                                readsbetw.add(s1);
                        }
                        if (readsbetw.size() == 2)
                        {
                            ArrayList<Read> pair = new ArrayList<Read>();
                            for (Read r : hap.reads)
                            {
                                StringTokenizer stk = new StringTokenizer(r.name,"_");
                                String stnm = stk.nextToken();
                                if (stnm.equalsIgnoreCase(readsbetw.get(0)))
                                {
                                    pair.add(r);
                                    break;
                                }
                            }
                            for (Read r : hap.reads)
                            {
                                StringTokenizer stk = new StringTokenizer(r.name,"_");
                                String stnm = stk.nextToken();
                                if (stnm.equalsIgnoreCase(readsbetw.get(1)))
                                {
                                    pair.add(r);
                                    break;
                                }
                            }
                            closepairs.add(pair);
                        }
                        if (readsbetw.size() > 2)
                        {
                            ArrayList<ArrayList<Read>> possibpairs = new ArrayList<ArrayList<Read>>();
                            for (int u = 0; u < readsbetw.size(); u++)
                                for (int v = u + 1; v < readsbetw.size(); v++)
                                {
                                    ArrayList<Read> pospair = new ArrayList<Read>();
                                    for (Read r : hap.reads)
                                    {
                                        StringTokenizer stk = new StringTokenizer(r.name,"_");
                                        String stnm = stk.nextToken();
                                        if (stnm.equalsIgnoreCase(readsbetw.get(u)))
                                        {
                                            pospair.add(r);
                                            break;
                                        }
                                    }
                                    for (Read r : hap.reads)
                                    {
                                        StringTokenizer stk = new StringTokenizer(r.name,"_");
                                        String stnm = stk.nextToken();
                                        if (stnm.equalsIgnoreCase(readsbetw.get(v)))
                                        {
                                            pospair.add(r);
                                            break;
                                        }
                                    }
                                    possibpairs.add(pospair);
                                }
                            ArrayList<Read> optpair = new ArrayList<Read>();
                            int optnucldiff = Integer.MAX_VALUE;
                            int opthomopdiff = Integer.MAX_VALUE;
                            double optdominparam = Double.MAX_VALUE;
                            for (ArrayList<Read> ar : possibpairs)
                            {
                                FileWriter fw_alin = new FileWriter("allign_input.fas");
                                fw_alin.write(">" + ar.get(0).name + "\n" + ar.get(0).nucl + "\n");
                                fw_alin.write(">" + ar.get(1).name + "\n" + ar.get(1).nucl + "\n");
                                fw_alin.close();
                                String param = "ClustalW2//clustalw2 -INFILE=" +"allign_input.fas" + " -OUTFILE=" + "allign_output.fas";
                                param+= " -OUTPUT=FASTA -DNAMATRIX=IUB -GAPOPEN=" + gapop +" -GAPEXT=" + gapext + " -TYPE=DNA -PWDNAMATRIX=IUB -PWGAPOPEN=" + gapop+ " -PWGAPEXT=" + gapext;
                                p = run.exec(param);
                                stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
                                while ((s = stdInput.readLine()) != null)
                                {
                                    System.out.println(s);
                                }
                                p.waitFor();
                                if (p.exitValue()!= 0)
                                {
                                    System.out.println("Error in allgnment program");
                                }
                                stdInput.close();
                                DataSet alignment = new DataSet("allign_output.fas",'c');
                                File f = new File("allign_input.fas");
                                f.delete();
                                f = new File("allign_output.fas");
                                f.delete();
                                f = new File("allign_input.dnd");
                                f.delete();
                                String fir = alignment.reads.get(0).getNucl();
                                String sec = alignment.reads.get(1).getNucl();
                                int startcompdiff = 0;
                                int endcompdiff = fir.length()-1;
                                String sz = "";
                                if (fir.charAt(0) == '-')
                                    sz = fir;
                                else
                                    sz = sec;
                                while (sz.charAt(startcompdiff) == '-')
                                    startcompdiff++;
                                if (fir.charAt(fir.length()-1) == '-')
                                    sz = fir;
                                else
                                    sz = sec;
                                while (sz.charAt(endcompdiff) == '-')
                                    endcompdiff--;
                                int diff = 0;
                                int homopdiff = 0;
                                for (i = startcompdiff; i <= endcompdiff; i++)
                                    if (fir.charAt(i)!= sec.charAt(i))
                                        diff++;
                                if (diff < optnucldiff)
                                {
                                    optnucldiff = diff;
                                    optpair = ar;
                                }
                            }
                            closepairs.add(optpair);
                        }
                    }
                }
                for (ArrayList<Read> ar : closepairs)
                {
                    System.out.println(ar.get(0).name);
                    System.out.println(ar.get(1).name);
                    System.out.println();
                }
                boolean correction = false;
                for (ArrayList<Read> ar : closepairs)
                {                    
                    Read large = null;
                    Read small = null;
                    if (ar.get(0).frequency >= ar.get(1).frequency)
                    {
                        large = ar.get(0);
                        small = ar.get(1);
                    }
                    else
                    {
                        large = ar.get(1);
                        small = ar.get(0);
                    }
//                    if (large.frequency < small.frequency * dominparam)
//                        continue;

                    run=Runtime.getRuntime();
                    p=null;
                    FileWriter fw_alin = new FileWriter("allign_input.fas");
                    fw_alin.write(">" + large.name + "\n" + large.nucl + "\n");
                    fw_alin.write(">" + small.name + "\n" + small.nucl + "\n");
                    fw_alin.close();
                    String param = "ClustalW2//clustalw2 -INFILE=" +"allign_input.fas" + " -OUTFILE=" + "allign_output.fas";
                    param+= " -OUTPUT=FASTA -DNAMATRIX=IUB -GAPOPEN=" + gapop +" -GAPEXT=" + gapext + " -TYPE=DNA -PWDNAMATRIX=IUB -PWGAPOPEN=" + gapop+ " -PWGAPEXT=" + gapext;
                    p = run.exec(param);
                    stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
                    while ((s = stdInput.readLine()) != null)
                    {
                        System.out.println(s);
                    }
                    p.waitFor();
                    if (p.exitValue()!= 0)
                    {
                        System.out.println("Error in allgnment program");
                    }
                    stdInput.close();
                    DataSet alignment = new DataSet("allign_output.fas",'c');
                    File f = new File("allign_input.fas");
                    f.delete();
                    f = new File("allign_output.fas");
                    f.delete();
                    f = new File("allign_input.dnd");
                    f.delete();
                    String large_align = "";
                    String small_align = "";
                    for (Read r: alignment.reads)
                    {
                        if (r.name.equalsIgnoreCase(large.name))
                            large_align = r.nucl;
                        if (r.name.equalsIgnoreCase(small.name))
                            small_align = r.nucl;
                    }
                    int l = large_align.length();
                    int large_pos = -1;
                    int small_pos = -1;

                    int startcompdiff = 0;
                    int endcompdiff = l-1;
                    String sz = "";
                    if (large_align.charAt(0) == '-')
                        sz = large_align;
                    else
                        sz = small_align;
                    while (sz.charAt(startcompdiff) == '-')
                        startcompdiff++;
                    if (large_align.charAt(l-1) == '-')
                        sz = large_align;
                    else
                        sz = small_align;
                    while (sz.charAt(endcompdiff) == '-')
                        endcompdiff--;

                    int diff = 0;
                    for (i = startcompdiff; i <= endcompdiff; i++)
                        if (large_align.charAt(i)!= small_align.charAt(i))
                            diff++;

                    if (diff <= nucldiffparam)
                    {
                        boolean corrhomop = false;
                        for (i = 0; i < startcompdiff; i++)
                        {
                            if (large_align.charAt(i)!='-')
                                large_pos++;
                            if (small_align.charAt(i)!='-')
                                small_pos++;
                        }
                        for (i = startcompdiff; i <= endcompdiff; i++)
                        {
                            if (large_align.charAt(i)!='-')
                                large_pos++;
                            if (small_align.charAt(i)!='-')
                                small_pos++;
                            if (small_pos == -1)
                                continue;
                            if ((large_align.charAt(i) == '-')&&(small_align.charAt(i)!= '-'))
                            {
                                if ((i > 0)&&(large_align.charAt(i-1) == small_align.charAt(i-1))&&(small_align.charAt(i-1) == small_align.charAt(i)))
                                {
                                    small.corrections[small_pos] = new Correction('D');
                                    corrhomop = true;
                                    correction = true;
                                }
                                if ((i < l-1)&&(large_align.charAt(i+1) == small_align.charAt(i+1))&&(small_align.charAt(i+1) == small_align.charAt(i)))
                                {
                                    small.corrections[small_pos] = new Correction('D');
                                    corrhomop = true;
                                    correction = true;
                                }
                            }
                            if ((large_align.charAt(i) != '-')&&(small_align.charAt(i)== '-'))
                            {
                                if ((i > 0)&&(large_align.charAt(i-1) == small_align.charAt(i-1))&&(large_align.charAt(i-1) == large_align.charAt(i)))
                                {
                                    small.corrections[small_pos] = new Correction('I',large_align.charAt(i));
                                    corrhomop = true;
                                    correction = true;
                                }
                                if ((i < l-1)&&(large_align.charAt(i+1) == small_align.charAt(i+1))&&(large_align.charAt(i+1) == large_align.charAt(i)))
                                {
                                    small.corrections[small_pos] = new Correction('I',large_align.charAt(i));
                                    corrhomop = true;
                                    correction = true;
                                }
                            }
                         }
                        if (corrhomop)
                            continue;
                        large_pos = -1;
                        small_pos = -1;
                        if (large.frequency >= small.frequency * dominparamonenucl)
                        {
                            for (i = 0; i < startcompdiff; i++)
                            {
                                if (large_align.charAt(i)!='-')
                                    large_pos++;
                                if (small_align.charAt(i)!='-')
                                    small_pos++;
                            }
                            for (i = startcompdiff; i <= endcompdiff; i++)
                            {
                                if (large_align.charAt(i)!='-')
                                    large_pos++;
                                if (small_align.charAt(i)!='-')
                                    small_pos++;
                                if (small_pos == -1)
                                    continue;
                                if ((large_align.charAt(i) == '-')&&(small_align.charAt(i)!= '-'))
                                {
                                    small.corrections[small_pos] = new Correction('D');
                                    correction = true;
                                }
                                if ((large_align.charAt(i) != '-')&&(small_align.charAt(i)== '-'))
                                {
                                    small.corrections[small_pos] = new Correction('I',large_align.charAt(i));
                                    correction = true;
                                }
                             }
                        }
                        continue;
                     }

                    if (large.frequency < small.frequency * dominparamgen)
                        continue;
                    boolean correctionthispair = false;
                    for (i = 0; i < l; i++)
                    {
                        if (large_align.charAt(i)!='-')
                            large_pos++;
                        if (small_align.charAt(i)!='-')
                            small_pos++;
                        if (small_pos == -1)
                            continue;
                        if ((large_align.charAt(i) == '-')&&(small_align.charAt(i)!= '-'))
                        {
                            if ((i > 0)&&(large_align.charAt(i-1) == small_align.charAt(i-1))&&(small_align.charAt(i-1) == small_align.charAt(i)))
                            {
                                small.corrections[small_pos] = new Correction('D');
                                correction = true;
                                correctionthispair = true;
                            }
                            if ((i < l-1)&&(large_align.charAt(i+1) == small_align.charAt(i+1))&&(small_align.charAt(i+1) == small_align.charAt(i)))
                            {
                                small.corrections[small_pos] = new Correction('D');
                                correction = true;
                                correctionthispair = true;
                            }
                        }
                        if ((large_align.charAt(i) != '-')&&(small_align.charAt(i)== '-'))
                        {
                            if ((i > 0)&&(large_align.charAt(i-1) == small_align.charAt(i-1))&&(large_align.charAt(i-1) == large_align.charAt(i)))
                            {
                                small.corrections[small_pos] = new Correction('I',large_align.charAt(i));
                                correction = true;
                                correctionthispair = true;
                            }
                            if ((i < l-1)&&(large_align.charAt(i+1) == small_align.charAt(i+1))&&(large_align.charAt(i+1) == large_align.charAt(i)))
                            {
                                small.corrections[small_pos] = new Correction('I',large_align.charAt(i));
                                correction = true;
                                correctionthispair = true;
                            }
                        }
                     }
                    if (correctionthispair)
                        continue;
                    int inddiff = 0;
                    for (i = startcompdiff; i <= endcompdiff; i++)
                        if (large_align.charAt(i)!= small_align.charAt(i))
                            if ((large_align.charAt(i)=='-')||(small_align.charAt(i)=='-'))
                                inddiff++;
                    if ((inddiff == diff)&&(large.frequency > small.frequency * dominparamgenins))
                    {
                        large_pos = -1;
                        small_pos = -1;
                        for (i = 0; i < l; i++)
                        {
                            if (large_align.charAt(i)!='-')
                                large_pos++;
                            if (small_align.charAt(i)!='-')
                                small_pos++;
                            if (small_pos == -1)
                                continue;
                            if ((large_align.charAt(i) == '-')&&(small_align.charAt(i)!= '-'))
                            {
                                    small.corrections[small_pos] = new Correction('D');
                                    correction = true;
                            }
                            if ((large_align.charAt(i) != '-')&&(small_align.charAt(i)== '-'))
                            {
                                    small.corrections[small_pos] = new Correction('I',large_align.charAt(i));
                                    correction = true;
                            }
                         }
                      }
                }
                DataSet corhap = new DataSet(hap);
                corhap.findHaplotypes();
                if (iter == niterpostp - 1)
                    addr = addr + "_PostprocPair.fas";
                StringTokenizer st1 = new StringTokenizer(addr,"_");
                String tag = st1.nextToken();
                corhap.PrintHaplotypesWithNameTag(addr, tag);
//                corhap.PrintHaplotypes(addr);
                System.out.println("");
                hap = new DataSet(addr,'c');
                //// ???????????????????????
                if (((hap.reads.size() == 1)&&(iter!= niterpostp - 1))|| (!correction))
                {
                   File f1 = new File(addr);
                   File f2 = new File(addr + "_PostprocPair.fas");
                   f1.renameTo(f2);
                   return;
                }
           }
/*           DataSet corhap = new DataSet(hap);
           corhap.findHaplotypes();
           corhap.PrintHaplotypes(addr+"_postprocessedPairwise.fas");
           System.out.println("Postprocessing finished!");*/
       }
       void postprocessFinal(String addr) throws IOException
       {
           File fl = new File(addr);
           if (!fl.exists())
           {
                System.out.println("No such file!");
                return;
           }
           if (ds.haplotypes.size() == 1)
           {
                System.out.println("No postprocessing needed");
                File f1 = new File(addr);
                File f2 = new File(addr + "_postprocessed.fas");
                f1.renameTo(f2);
                return;
           }
           DataSet dats = new DataSet(addr);
           if (dats.getNreads() == 1)
           {
                System.out.println("No postprocessing needed");
                return;
           }
           StrComp sc = new StrComp();
           Collections.sort(dats.reads, sc);
           int m = dats.reads.size()/2;
           int median = dats.reads.get(m).getLength();
           System.out.println((double) dats.reads.get(dats.reads.size()-1).getLength() / median);
           for (Read r : dats.reads)
               System.out.println(r.getLength());
           System.out.println();


       }
       void calcStatHomop (String addr) throws IOException
       {
           FileReader fr = new FileReader(addr);
           BufferedReader br = new BufferedReader(fr);
           String name = br.readLine();
           String s = br.readLine();
           FileWriter fw = new FileWriter(addr+"_HomopStat.txt");
           fw.write("name ");
           for (int i = 1; i <= 10; i++)
               fw.write("homop" + i + " ");
           fw.write("\n");
           while(s!=null)
           {
               int[] len = new int[11];
               int i = 0;
               while (i < s.length())
               {
                   int l = 0;
                   int j = i;
                   while ((j < s.length())&&(s.charAt(j)==s.charAt(i)))
                   {
                           j++;
                           l++;
                   }
                   i = j;
                   len[l]++;
               }
               fw.write(name + " ");
               for (i = 1; i <= 10; i++)
                   fw.write(len[i] + " ");
               fw.write("\n");
               name = br.readLine();
               s = br.readLine();
           }
           fw.close();
           fr.close();

       }
       void postprocessHaplotypesManualAllignedET(String addr, ArrayList<Integer> positions,ArrayList<Character> type, ArrayList<Character> nucleotide) throws IOException, InterruptedException
        {
                DataSet hap = new DataSet(addr,'c');
                if (hap.getNreads() == 1)
                {
                    System.out.println("No postprocessing needed");
                    return;
                }
                
                for (Read r : hap.reads)
                {
                    int ndeleted = 0;
                    for (int i = 0; i < positions.size(); i++)
                    {
                        if (type.get(i) == 'D')
                        {
                            r.removeNucl(positions.get(i) - 1 - ndeleted);
                            ndeleted++;
                        }
                        if (type.get(i) == 'I')
                            if (r.getNucl((positions.get(i) - 1 - ndeleted)) == '-')
                                r.replaceNucl(positions.get(i) - 1 - ndeleted, nucleotide.get(i));
                    }
                }
                
                for (Read r: hap.reads)
                {
                    StringTokenizer st = new StringTokenizer(r.name,"_");
                    st.nextToken();
                    st.nextToken();
                    st.nextToken();
                    st.nextToken();
                    st.nextToken();
                    r.frequency = Integer.parseInt(st.nextToken());
                }
                
              DataSet corhap = new DataSet(hap);
              corhap.findHaplotypes();
              StringTokenizer st1 = new StringTokenizer(addr,"_");
              String tag = st1.nextToken();
              corhap.PrintHaplotypesWithNameTag(addr+"_postprocessed.fas", tag);
        }
       void postprocessHaplotypesNuclSwitch(String addr, String ref_addr, String idmethod, double gapop, double gapext) throws FileNotFoundException, IOException
       // <editor-fold defaultstate="collapsed" desc=" DESCRIPTION ">
       {
           DataSet hapl = new DataSet(addr,"ET");
           Read ref = (new DataSet(ref_addr)).reads.get(0);
           int count = 0;
           for (Read h : hapl.reads)
           {
                    count++;
                    System.out.println("Aligning " + count + File.separator + hapl.reads.size());
                    Runtime run=Runtime.getRuntime();
                    Process p=null;
                    FileWriter fw_alin = new FileWriter("allign_input.fas");
                    fw_alin.write(">" + "reference" + "\n" + ref.nucl + "\n");
                    fw_alin.write(">" + h.name + "\n" + h.nucl + "\n");
                    fw_alin.close();
                    String param = "ClustalW2//clustalw2 -INFILE=" +"allign_input.fas" + " -OUTFILE=" + "allign_output.fas";
                    param+= " -OUTPUT=FASTA -DNAMATRIX=IUB -GAPOPEN=" + gapop +" -GAPEXT=" + gapext + " -TYPE=DNA -PWDNAMATRIX=IUB -PWGAPOPEN=" + gapop+ " -PWGAPEXT=" + gapext;
                    p = run.exec(param);
                    BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
                    String s;
                    while ((s = stdInput.readLine()) != null)
                    {
//                        System.out.println(s);
                    }
                    try {
                        p.waitFor();
                    } catch (InterruptedException ex) {
                         //Logger.getLogger(Corrector.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    if (p.exitValue()!= 0)
                    {
                        System.out.println("Error in allgnment program");
                    }
                    stdInput.close();
                    DataSet alignment = new DataSet("allign_output.fas",'c');
                    File f = new File("allign_input.fas");
                    f.delete();
                    f = new File("allign_output.fas");
                    f.delete();
                    f = new File("allign_input.dnd");
                    f.delete();
                    String ref_align = "";
                    String hapl_align = "";
                    for (Read r: alignment.reads)
                    {
                        if (r.name.equalsIgnoreCase("reference"))
                            ref_align = r.nucl;
                        if (r.name.equalsIgnoreCase(h.name))
                            hapl_align = r.nucl;
                    }
                    int l = ref_align.length();
                    int hapl_pos = -1;
                    for (int i = 0; i < l; i++)
                    {
                        if (hapl_align.charAt(i)!='-')
                                hapl_pos++;
                        if ((ref_align.charAt(i) != hapl_align.charAt(i)) && (ref_align.charAt(i) != '-') && (hapl_align.charAt(i) != '-'))
                        {
                            // looking forward
                            int j = i;
                            int hapl_pos1 = hapl_pos;
                            while ((ref_align.charAt(j) != hapl_align.charAt(i)) && (j < l))
                            {
                                j++;
                                if (hapl_align.charAt(j) != '-')
                                    hapl_pos1++;
                            }
                            if (j == l)
                                continue;
                            String cand = hapl_align.substring(i+1, j+1) + hapl_align.charAt(i);
                            if (cand.equalsIgnoreCase(ref_align.substring(i, j+1)))
                            {
                                h.corrections[hapl_pos] = new Correction('D');
                                h.corrections[hapl_pos1] = new Correction('I',hapl_align.charAt(i));
                                continue;
                            }
                            
                            // looking backward
                            
                            j = i;
                            hapl_pos1 = hapl_pos;
                            while ((ref_align.charAt(j) != hapl_align.charAt(i)) && (j >= 0))
                            {
                                j--;
                                if (hapl_align.charAt(j) != '-')
                                    hapl_pos1--;
                            }
                            if ((j == -1) || (hapl_pos1 == 0))
                                continue;
                            cand = "" + hapl_align.charAt(i) + hapl_align.substring(j, i);
                            if (cand.equalsIgnoreCase(ref_align.substring(j, i+1)))
                            {
                                h.corrections[hapl_pos] = new Correction('D');
                                h.corrections[hapl_pos1-1] = new Correction('I',hapl_align.charAt(i));
                                continue;
                            }
                                                        
                        }
                        
                    }
           }
           DataSet corhapl = new DataSet(hapl);
           corhapl.findHaplotypes();
           corhapl.PrintHaplotypes(addr + "_PostprocShift.fas");
       }
       // </editor-fold>
       void postprocessHaplotypesNuclSwitch(String addr, String idmethod, double gapop, double gapext, int homop_size) throws FileNotFoundException, IOException
       // <editor-fold defaultstate="collapsed" desc=" DESCRIPTION ">
       {
           DataSet hapl = null;
           if (idmethod.equalsIgnoreCase("RAW"))
               hapl = new DataSet(addr);
           else
               hapl = new DataSet(addr,idmethod);
           Read ref = hapl.findMostFreqRead();
           int count = 0;
           for (Read h : hapl.reads)
           {
                    count++;
                    System.out.println("Aligning " + count + "\\" + hapl.reads.size());
                    Runtime run=Runtime.getRuntime();
                    Process p=null;
                    FileWriter fw_alin = new FileWriter("allign_input.fas");
                    fw_alin.write(">" + "reference" + "\n" + ref.nucl + "\n");
                    fw_alin.write(">" + h.name + "\n" + h.nucl + "\n");
                    fw_alin.close();
                    String param = "ClustalW2//clustalw2 -INFILE=" +"allign_input.fas" + " -OUTFILE=" + "allign_output.fas";
                    param+= " -OUTPUT=FASTA -DNAMATRIX=IUB -GAPOPEN=" + gapop +" -GAPEXT=" + gapext + " -TYPE=DNA -PWDNAMATRIX=IUB -PWGAPOPEN=" + gapop+ " -PWGAPEXT=" + gapext;
                    p = run.exec(param);
                    BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
                    String s;
                    while ((s = stdInput.readLine()) != null)
                    {
//                        System.out.println(s);
                    }
            try {
                p.waitFor();
            } catch (InterruptedException ex) {
               // Logger.getLogger(Corrector.class.getName()).log(Level.SEVERE, null, ex);
            }
                    if (p.exitValue()!= 0)
                    {
                        System.out.println("Error in allgnment program");
                    }
                    stdInput.close();
                    DataSet alignment = new DataSet("allign_output.fas",'c');
                    File f = new File("allign_input.fas");
                    f.delete();
                    f = new File("allign_output.fas");
                    f.delete();
                    f = new File("allign_input.dnd");
                    f.delete();
                    String ref_align = "";
                    String hapl_align = "";
                    for (Read r: alignment.reads)
                    {
                        if (r.name.equalsIgnoreCase("reference"))
                            ref_align = r.nucl;
                        if (r.name.equalsIgnoreCase(h.name))
                            hapl_align = r.nucl;
                    }
                    int l = ref_align.length();
                    int hapl_pos = -1;
                    for (int i = 0; i < l; i++)
                    {
                        if (hapl_align.charAt(i)!='-')
                                hapl_pos++;
                        if ((i == 0) || (i == l - 1))
                            continue;
                        
                        if ((ref_align.charAt(i) != hapl_align.charAt(i)) && (ref_align.charAt(i) != '-') && (hapl_align.charAt(i) != '-'))
                        {
                            // looking forward
                            int j = i;
                            int hapl_pos1 = hapl_pos;
                            while ((j < l) && (ref_align.charAt(j) != hapl_align.charAt(i)))
                            {
                                j++;
                                if (j == l)
                                    break;
                                if (hapl_align.charAt(j) != '-')
                                    hapl_pos1++;
                            }
                            if (j == l)
                                continue;
                            String cand = hapl_align.substring(i+1, j+1) + hapl_align.charAt(i);
                            if (cand.equalsIgnoreCase(ref_align.substring(i, j+1)))
                            {
                                boolean wasCorr = false;
                                for (int u = hapl_pos; u <= hapl_pos1; u++)
                                    if (h.corrections[u] != null)
                                    {
                                        wasCorr = true;
                                        break;
                                    }
                                if (!wasCorr)
                                {
                                    h.corrections[hapl_pos] = new Correction('D');
                                    h.corrections[hapl_pos1] = new Correction('I',hapl_align.charAt(i));
                                    continue;
                                }
                            }
                            
                            // looking backward
                            
                            j = i;
                            hapl_pos1 = hapl_pos;
                            while ((j >= 0) && (ref_align.charAt(j) != hapl_align.charAt(i)))
                            {
                                j--;
                                if (j < 0)
                                    break;
                                if (hapl_align.charAt(j) != '-')
                                    hapl_pos1--;
                            }
                            if ((j == -1) || (hapl_pos1 == 0))
                                continue;
                            cand = "" + hapl_align.charAt(i) + hapl_align.substring(j, i);
                            if (cand.equalsIgnoreCase(ref_align.substring(j, i+1)))
                            {
                                boolean wasCorr = false;
                                for (int u = hapl_pos1 - 1; u <= hapl_pos; u++)
                                    if (h.corrections[u] != null)
                                    {
                                        wasCorr = true;
                                        break;
                                    }
                                if (!wasCorr)
                                {
                                    h.corrections[hapl_pos] = new Correction('D');
                                    h.corrections[hapl_pos1-1] = new Correction('I',hapl_align.charAt(i));
                                    continue;
                                }
                            }
                                                        
                        }
                        // homopolymer-related swithches
                        
                        if ((ref_align.charAt(i) == '-') && (hapl_align.charAt(i) != '-'))
                        {
                            char susp_nucl = hapl_align.charAt(i);
                            
                            // looking backward
                            
                            int j = i - 1;
                            while ((j > 1) && (ref_align.charAt(j) != susp_nucl) && (hapl_align.charAt(j) != susp_nucl))
                                j--;                            
                            if (j <= 1)
                                continue;
                            if (j == i-1)
                                continue;
                            boolean inCycle = true;
                            if (!ref_align.substring(j+1, i).equalsIgnoreCase(hapl_align.substring(j+1, i)))
                                inCycle = false;
                            HashSet<Character> hs = new HashSet<Character>();
                            String contract_homop ="" + ref_align.charAt(j+1);
                            for (int u = j+2; u <= i-1; u++)
                                if (ref_align.charAt(u) != contract_homop.charAt(contract_homop.length() - 1))
                                    contract_homop = contract_homop + ref_align.charAt(u);
                            for (int u = 0; u < contract_homop.length(); u++)
                                if (hs.contains(contract_homop.charAt(u)))
                                {
                                    inCycle = false;
                                    break;
                                }
                                else
                                    hs.add(contract_homop.charAt(u));
                            if (inCycle)
                            {
                                hs = new HashSet<Character>();
                                hs.add(susp_nucl);
                                hs.add('-');
                                int u = j;
                                while ((u >= 0) && (hs.contains(ref_align.charAt(u))) && (hs.contains(hapl_align.charAt(u))))
                                    u--;
                                if (j - u >= homop_size)
                                {
                                    h.corrections[hapl_pos] = new Correction('D');
                                    continue;
                                }
                            }
                            
                             // looking forward
                            
                            j = i + 1;
                            while ((j <= l-2) && (ref_align.charAt(j) != susp_nucl) && (hapl_align.charAt(j) != susp_nucl))
                                j++;                            
                            if (j >= l-1)
                                continue;
                            if (j == i+1)
                                continue;
                            inCycle = true;
                            if (!ref_align.substring(i+1, j).equalsIgnoreCase(hapl_align.substring(i+1, j)))
                                inCycle = false;
                            hs = new HashSet<Character>();    
                            contract_homop ="" + ref_align.charAt(i+1);
                            for (int u = i+2; u <= j-1; u++)
                                if (ref_align.charAt(u) != contract_homop.charAt(contract_homop.length() - 1))
                                    contract_homop = contract_homop + ref_align.charAt(u);
                            for (int u = 0; u < contract_homop.length(); u++)
                                if (hs.contains(contract_homop.charAt(u)))
                                {
                                    inCycle = false;
                                    break;
                                }
                                else
                                    hs.add(contract_homop.charAt(u));
                            if (inCycle)
                            {                                 
                                hs = new HashSet<Character>();
                                hs.add(susp_nucl);
                                hs.add('-');
                                int u = j;
                                while ((u < l) && (hs.contains(ref_align.charAt(u))) && (hs.contains(hapl_align.charAt(u))))
                                    u++;
                                if (u - j >= homop_size)
                                    h.corrections[hapl_pos] = new Correction('D');
                            }
                            
                        }
                        
                    }
           }
           DataSet corhapl = new DataSet(hapl);
//           corhapl.PrintCorrectedReads(addr + "_correctedShiftReads.fas");
           corhapl.findHaplotypes();
           corhapl.PrintHaplotypes(addr + "_PostprocShift.fas");
       }
       // </editor-fold>
       void delGaps(String addr) throws FileNotFoundException, IOException
    {
                FileReader fr = new FileReader(addr);
                FileWriter fw = new FileWriter(addr +"_delGaps.fas");
                BufferedReader br = new BufferedReader(fr);
                String s = br.readLine();
                while(s.equalsIgnoreCase(""))
                    s = br.readLine();
		String nucl = "";
		String name = s;
		s = br.readLine();
		while (s!= null)
		{
			if (s.length() == 0)
			{
				s = br.readLine();
				continue;
			}
			if(s.charAt(0) == '>')
			{
                                        fw.write(name + "\n");
                                        String nucl1 = "";
                                        for (int i = 0; i < nucl.length(); i++)
                                            if (nucl.charAt(i)!='-')
                                                nucl1 = nucl1 + nucl.charAt(i);
                                        fw.write(nucl1 + "\n");
					nucl = "";
					name = s;
			}
			else
				nucl+=s.toUpperCase();
			s = br.readLine();
		}
                fw.write(name + "\n");
                String nucl1 = "";
                for (int i = 0; i < nucl.length(); i++)
                    if (nucl.charAt(i)!='-')
                        nucl1 = nucl1 + nucl.charAt(i);
                fw.write(nucl1 + "\n");
                br.close();
                fr.close();
                fw.close();
        }
       int findClosest(String refer, double gapop, double gapext) throws IOException
       {
           int min = 100000;
           String minName = "";
           int count = 0;
           for (Read r : ds.reads)
           {
                    count++;
                    System.out.println("Aligning " + count + "\\" + ds.reads.size());
                    if (r.getLength() < 200)
                        continue;
                    Runtime run=Runtime.getRuntime();
                    Process p=null;
                    FileWriter fw_alin = new FileWriter("allign_input.fas");
                    fw_alin.write(">" + "reference" + "\n" + refer + "\n");
                    fw_alin.write(">" + r.name + "\n" + r.nucl + "\n");
                    fw_alin.close();
                    String param = "ClustalW2//clustalw2 -INFILE=" +"allign_input.fas" + " -OUTFILE=" + "allign_output.fas";
                    param+= " -OUTPUT=FASTA -DNAMATRIX=IUB -GAPOPEN=" + gapop +" -GAPEXT=" + gapext + " -TYPE=DNA -PWDNAMATRIX=IUB -PWGAPOPEN=" + gapop+ " -PWGAPEXT=" + gapext;
                    p = run.exec(param);
                    BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
                    String s;
                    while ((s = stdInput.readLine()) != null)
                    {
//                        System.out.println(s);
                    }
                    if (p.exitValue()!= 0)
                    {
                        System.out.println("Error in allgnment program");
                    }
                    DataSet alignment = new DataSet("allign_output.fas",'c');
                    File f = new File("allign_input.fas");
                    f.delete();
                    f = new File("allign_output.fas");
                    f.delete();
                    f = new File("allign_input.dnd");
                    f.delete();
                    String ref_align = "";
                    String read_align = "";
                    for (Read r1: alignment.reads)
                    {
                        if (r1.name.equalsIgnoreCase("reference"))
                            ref_align = r1.nucl;
                        if (r1.name.equalsIgnoreCase(r.name))
                            read_align = r1.nucl;
                    }
                    int start = 0;
                    int end = ref_align.length() - 1;
                    while ((ref_align.charAt(start)=='-') || (read_align.charAt(start))=='-')
                        start++;
                    while ((ref_align.charAt(end)=='-') || (read_align.charAt(end))=='-')
                        end--;
                    int dist = 0;
                    for (int i = start; i <= end; i++)
                        if (ref_align.charAt(i) != read_align.charAt(i))
                            dist++;
                    if (dist < min)
                    {
                        min = dist;
                        minName = r.name;
                    }
           }
           System.out.println(min);
           System.out.println(minName);
           return min;
       }
        int findClosestHapl(String refer, double gapop, double gapext) throws IOException
                // <editor-fold defaultstate="collapsed" desc=" DESCRIPTION ">
       {
           int min = 100000;
           String minName = "";
           int count = 0;
           for (Haplotype r : ds.haplotypes)
           {
                    count++;
                    System.out.println("Aligning " + count + "\\" + ds.haplotypes.size());
                    Runtime run=Runtime.getRuntime();
                    Process p=null;
                    FileWriter fw_alin = new FileWriter("allign_input.fas");
                    fw_alin.write(">" + "reference" + "\n" + refer + "\n");
                    fw_alin.write(">" + r.name + "\n" + r.nucl + "\n");
                    fw_alin.close();
                    String param = "ClustalW2//clustalw2 -INFILE=" +"allign_input.fas" + " -OUTFILE=" + "allign_output.fas";
                    param+= " -OUTPUT=FASTA -DNAMATRIX=IUB -GAPOPEN=" + gapop +" -GAPEXT=" + gapext + " -TYPE=DNA -PWDNAMATRIX=IUB -PWGAPOPEN=" + gapop+ " -PWGAPEXT=" + gapext;
                    p = run.exec(param);
                    BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
                    String s;
                    while ((s = stdInput.readLine()) != null)
                    {
//                        System.out.println(s);
                    }
                    if (p.exitValue()!= 0)
                    {
                        System.out.println("Error in allgnment program");
                    }
                    DataSet alignment = new DataSet("allign_output.fas",'c');
                    File f = new File("allign_input.fas");
                    f.delete();
                    f = new File("allign_output.fas");
                    f.delete();
                    f = new File("allign_input.dnd");
                    f.delete();
                    String ref_align = "";
                    String read_align = "";
                    for (Read r1: alignment.reads)
                    {
                        if (r1.name.equalsIgnoreCase("reference"))
                            ref_align = r1.nucl;
                        if (r1.name.equalsIgnoreCase(r.name))
                            read_align = r1.nucl;
                    }
                    int start = 0;
                    int end = ref_align.length() - 1;
                    while ((ref_align.charAt(start)=='-') || (read_align.charAt(start))=='-')
                        start++;
                    while ((ref_align.charAt(end)=='-') || (read_align.charAt(end))=='-')
                        end--;
                    int dist = 0;
                    for (int i = start; i <= end; i++)
                        if (ref_align.charAt(i) != read_align.charAt(i))
                            dist++;
                    if (dist < min)
                    {
                        min = dist;
                        minName = r.name;
                    }
           }
           System.out.println(min);
           return min;
       }
               // </editor-fold>
}
