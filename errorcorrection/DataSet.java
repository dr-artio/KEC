package ErrorCorrection;
import java.util.*;
import java.io.*;
import java.math.BigDecimal;

import static ErrorCorrection.Poisson.*;
import static ErrorCorrection.Exponent.*;
import java.util.logging.Level;
import java.util.logging.Logger;

public class DataSet {
	int freqThr;
	int lengthThr;
	int additionalfreqThrMult;
	int maxAllErrorsPerc;
	public ArrayList<Read> reads;
	ArrayList<Read> badReads;
	ArrayList<Read> goodReads;
	ArrayList<Kmer_general> allKmers;
	ArrayList<Haplotype> haplotypes;
	int k;
	String file_name;
	int[] valuesDistribution;
	int[] valuesDistribution1;
        int finderrorsseglen;
	public DataSet (String addr) throws IOException
	{
		lengthThr = 50;
		additionalfreqThrMult = 10;
		maxAllErrorsPerc = 40;
                freqThr = -1;
		reads = new ArrayList<Read>();
		allKmers = new ArrayList<Kmer_general>();
		haplotypes = new ArrayList<Haplotype>();
		file_name = addr;
                finderrorsseglen = k;
		FileReader fr = new FileReader(addr);
		BufferedReader br = new BufferedReader(fr);
		
		String s = br.readLine();
                while(s.equalsIgnoreCase(""))
                    s = br.readLine();

		String nucl = new String();
		String name = new String();
		name = "";
		int i = 1;
		int j = s.length();
		while ((i < s.length())&&(s.charAt(i)!=' '))
		{
			name+=s.charAt(i);
			i++;
		}
		s = br.readLine();
		
		HashMap<String,Integer> seq = new HashMap<String, Integer>();
                
		
		try
		{
		while (s!= null)
		{
//			System.out.println("Reading " + name);
			
			if (s.length() == 0)
			{
				s = br.readLine();
				continue;
			}
			if(s.charAt(0) == '>')
			{
				//	reads.add(new Read(nucl, name));
				
//					System.out.println(name);
				
					if (seq.containsKey(nucl))
						seq.put(nucl, seq.get(nucl)+1);
					else
						seq.put(nucl, 1);
                                        
					nucl = "";
					name = "";
					i = 1;
					while ((i < s.length())&&(s.charAt(i)!=' '))
					{
						name+=s.charAt(i);
						i++;
					}
			}
			else
				nucl+=s.toUpperCase();
			s = br.readLine();
		}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	//	reads.add(new Read(nucl, name));
		
		if (seq.containsKey(nucl))
			seq.put(nucl, seq.get(nucl)+1);
		else
			seq.put(nucl, 1);		
		Iterator it = seq.entrySet().iterator();
		int count = 1;
		while (it.hasNext())
		{
			Map.Entry me = (Map.Entry)it.next();
			String nucleo = (String) me.getKey(); 
			int f = (Integer) me.getValue();
			String nm = "read" + count;
			if (nucleo.length() > lengthThr)
				reads.add(new Read(nucleo,nm,f));
			count++;
		}
		
		fr.close();
	}
        public DataSet (String addr, int lt) throws IOException
	{
		lengthThr = lt; // 50
		additionalfreqThrMult = 10;
		maxAllErrorsPerc = 40;
                freqThr = -1;
		reads = new ArrayList<Read>();
		allKmers = new ArrayList<Kmer_general>();
		haplotypes = new ArrayList<Haplotype>();
		file_name = addr;
                finderrorsseglen = k;
		FileReader fr = new FileReader(addr);
		BufferedReader br = new BufferedReader(fr);

                String s = br.readLine();
                while(s.equalsIgnoreCase(""))
                    s = br.readLine();


		String nucl = new String();
		String name = new String();
		name = "";
		int i = 1;
		int j = s.length();
		while ((i < s.length())&&(s.charAt(i)!=' '))
		{
			name+=s.charAt(i);
			i++;
		}
		s = br.readLine();

		HashMap<String,Integer> seq = new HashMap<String, Integer>();

		try
		{
		while (s!= null)
		{
//			System.out.println("Reading " + name);

			if (s.length() == 0)
			{
				s = br.readLine();
				continue;
			}
			if(s.charAt(0) == '>')
			{
				//	reads.add(new Read(nucl, name));

//					System.out.println(name);

					if (seq.containsKey(nucl))
						seq.put(nucl, seq.get(nucl)+1);
					else
						seq.put(nucl, 1);

					nucl = "";
					name = "";
					i = 1;
					while ((i < s.length())&&(s.charAt(i)!=' '))
					{
						name+=s.charAt(i);
						i++;
					}
			}
			else
				nucl+=s.toUpperCase();
			s = br.readLine();
		}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	//	reads.add(new Read(nucl, name));

		if (seq.containsKey(nucl))
			seq.put(nucl, seq.get(nucl)+1);
		else
			seq.put(nucl, 1);
		Iterator it = seq.entrySet().iterator();
		int count = 1;
		while (it.hasNext())
		{
			Map.Entry me = (Map.Entry)it.next();
			String nucleo = (String) me.getKey();
			int f = (Integer) me.getValue();
			String nm = "read" + count;
			if (nucleo.length() > lengthThr)
				reads.add(new Read(nucleo,nm,f));
			count++;
		}

		fr.close();
	}
        public DataSet (String addr, int minlen, int maxlen) throws IOException
	{
		additionalfreqThrMult = 10;
		maxAllErrorsPerc = 40;
                freqThr = -1;
		reads = new ArrayList<Read>();
		allKmers = new ArrayList<Kmer_general>();
		haplotypes = new ArrayList<Haplotype>();
		file_name = addr;
                finderrorsseglen = k;
		FileReader fr = new FileReader(addr);
		BufferedReader br = new BufferedReader(fr);
		
		String s = br.readLine();
                while(s.equalsIgnoreCase(""))
                    s = br.readLine();

		String nucl = new String();
		String name = new String();
		name = "";
		int i = 1;
		int j = s.length();
		while ((i < s.length())&&(s.charAt(i)!=' '))
		{
			name+=s.charAt(i);
			i++;
		}
		s = br.readLine();
		
		HashMap<String,Integer> seq = new HashMap<String, Integer>();
                
		
		try
		{
		while (s!= null)
		{
//			System.out.println("Reading " + name);
			
			if (s.length() == 0)
			{
				s = br.readLine();
				continue;
			}
			if(s.charAt(0) == '>')
			{
				//	reads.add(new Read(nucl, name));
				
//					System.out.println(name);
				
					if (seq.containsKey(nucl))
						seq.put(nucl, seq.get(nucl)+1);
					else
						seq.put(nucl, 1);
                                        
					nucl = "";
					name = "";
					i = 1;
					while ((i < s.length())&&(s.charAt(i)!=' '))
					{
						name+=s.charAt(i);
						i++;
					}
			}
			else
				nucl+=s.toUpperCase();
			s = br.readLine();
		}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	//	reads.add(new Read(nucl, name));
		
		if (seq.containsKey(nucl))
			seq.put(nucl, seq.get(nucl)+1);
		else
			seq.put(nucl, 1);		
		Iterator it = seq.entrySet().iterator();
		int count = 1;
		while (it.hasNext())
		{
			Map.Entry me = (Map.Entry)it.next();
			String nucleo = (String) me.getKey(); 
			int f = (Integer) me.getValue();
			String nm = "read" + count;
			if ((nucleo.length() >= minlen)&&(nucleo.length() <= maxlen))
				reads.add(new Read(nucleo,nm,f));
			count++;
		}
		
		fr.close();
	}
	public DataSet (String addr, char c) throws IOException
	{
            // keep read name
		lengthThr = 50;
		additionalfreqThrMult = 3;
		maxAllErrorsPerc = 40;
                finderrorsseglen = k;
		reads = new ArrayList<Read>();
		allKmers = new ArrayList<Kmer_general>();
		haplotypes = new ArrayList<Haplotype>();
		file_name = addr;
		FileReader fr = new FileReader(addr);
		BufferedReader br = new BufferedReader(fr);
		
		String s = br.readLine();
                while(s.equalsIgnoreCase(""))
                    s = br.readLine();
		String nucl = new String();
		String name = new String();
		name = "";
		int i = 1;
		int j = s.length();
		while ((i < s.length())&&(s.charAt(i)!=' '))
		{
			name+=s.charAt(i);
			i++;
		}
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
					reads.add(new Read(nucl, name));	
					nucl = "";
					name = "";
					i = 1;
					while ((i < s.length())&&(s.charAt(i)!=' '))
					{
						name+=s.charAt(i);
						i++;
					}
			}
			else
				nucl+=s.toUpperCase();
			s = br.readLine();
		}
		reads.add(new Read(nucl, name));
		fr.close();
	}
	public DataSet(DataSet ds)
	{
		freqThr = ds.freqThr;
		lengthThr = ds.lengthThr;
		additionalfreqThrMult = ds.additionalfreqThrMult;
		maxAllErrorsPerc = ds.maxAllErrorsPerc;
		reads = new ArrayList<Read>();
		allKmers = new ArrayList<Kmer_general>();
		haplotypes = new ArrayList<Haplotype>();
		file_name = ds.file_name;
		k = ds.k;
                finderrorsseglen =ds.finderrorsseglen;
		
		HashMap<String,Integer> seq = new HashMap<String, Integer>();
		
		for(Read r : ds.reads)
		{
			String nucl = "";
			for (int i = 0; i < r.getLength(); i++)
				if (r.corrections[i] != null)
				{
					if (r.corrections[i].type == 'R')
						nucl = nucl + r.corrections[i].replacement;
					if (r.corrections[i].type == 'I')
						nucl = nucl + r.nucl.charAt(i) + r.corrections[i].replacement;
				}
				else
					nucl = nucl + r.nucl.charAt(i);
                        String nucl1 = "";
                        for (int i = 0; i < nucl.length(); i++)
                            if (nucl.charAt(i)!='-')
                                nucl1 = nucl1 + nucl.charAt(i);
                        nucl = nucl1;
			if (nucl.length() > lengthThr)
				if (seq.containsKey(nucl))
					seq.put(nucl, seq.get(nucl)+r.frequency);
				else
					seq.put(nucl, r.frequency);	
		}
		Iterator it = seq.entrySet().iterator();
		int count = 1;
		while (it.hasNext())
		{
			Map.Entry me = (Map.Entry)it.next();
			String nucleo = (String) me.getKey(); 
			int f = (Integer) me.getValue();
			String nm = "read" + count;
			if (nucleo.length() > lengthThr)
				reads.add(new Read(nucleo,nm,f));
			count++;
		}
	}
        public DataSet(String addr, String idmethod) throws FileNotFoundException, IOException
        {
                reads = new ArrayList<Read>();
		haplotypes = new ArrayList<Haplotype>();
                
                FileReader fr = new FileReader(addr);
		BufferedReader br = new BufferedReader(fr);
		
		String s = br.readLine();
                while(s.equalsIgnoreCase(""))
                    s = br.readLine();
		String nucl = new String();
		String name = new String();
		name = "";
		int i = 1;
		int j = s.length();
		while ((i < s.length())&&(s.charAt(i)!=' '))
		{
			name+=s.charAt(i);
			i++;
		}
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
                                        int f = 0;
                                        if (idmethod.equalsIgnoreCase("ET"))
                                        {
                                            StringTokenizer st = new StringTokenizer(name,"_");
                                            int nnn = st.countTokens();
                                            for (int tok = 1; tok <= nnn-1; tok++)
                                                st.nextToken();
                                            f = Integer.parseInt(st.nextToken());
                                        }
                                        
                                        if (idmethod.equalsIgnoreCase("KEC"))
                                        {
                                            StringTokenizer st = new StringTokenizer(name,"|");
                                            st.nextToken();
                                            f = (int) (10000*Double.parseDouble(st.nextToken()));
                                        }
                                                    
					reads.add(new Read(nucl, name,f));	
					nucl = "";
					name = "";
					i = 1;
					while ((i < s.length())&&(s.charAt(i)!=' '))
					{
						name+=s.charAt(i);
						i++;
					}
			}
			else
				nucl+=s.toUpperCase();
			s = br.readLine();
		}
                int f = 0;
                if (idmethod.equalsIgnoreCase("ET"))
                {
                    StringTokenizer st = new StringTokenizer(name,"_");
                    int nnn = st.countTokens();
                    for (int tok = 1; tok <= nnn-1; tok++)
                        st.nextToken();
                    f = Integer.parseInt(st.nextToken());
                }                
                if (idmethod.equalsIgnoreCase("KEC"))
                {
                    StringTokenizer st = new StringTokenizer(name,"|");
                    st.nextToken();
                    f = (int) (10000*Double.parseDouble(st.nextToken()));
                }
                
		reads.add(new Read(nucl, name,f));
		fr.close();           
        }
	public void setK(int i)
	{
		k=i;
	}
	void setFreqThr(int i)
	{
		freqThr = i;
	}
	public void setLenThr(int i)
	{
		lengthThr = i;
	}
	void setMaxAllErrorsPerc(int i)
	{
		maxAllErrorsPerc = i;
	}
	public void calculateKMersAndKCounts()
	{
		HashMap<String, ArrayList<Kmer_occur>> kmers= new HashMap<String, ArrayList<Kmer_occur>>();
		for (int j = 0; j < reads.size(); j++)
		{
			for (int i = 0; i< reads.get(j).nucl.length()-k+1; i++)
			{
				String s = reads.get(j).nucl.substring(i, i+k);
				if (kmers.containsKey(s))
					((ArrayList<Kmer_occur>) kmers.get(s)).add(new Kmer_occur(reads.get(j),j,i,k));
				else
				{
					ArrayList<Kmer_occur> a = new ArrayList<Kmer_occur>();
					a.add(new Kmer_occur(reads.get(j),j,i,k));
					kmers.put(s, a);
				}
			}
		}
		Iterator it = kmers.entrySet().iterator();
		int maxKCount = 0;
		while (it.hasNext())
		{
			Map.Entry me = (Map.Entry)it.next();
			String nucl =(String) me.getKey();
			ArrayList<Kmer_occur> readsOfKmer = ((ArrayList<Kmer_occur>) me.getValue());
			int kCount = 0;
			int kCount1 = 0;
//			HashSet<String> readsWhichContKmer = new HashSet<String>();
			for (int i = 0; i < readsOfKmer.size(); i++)
			{
				kCount+= readsOfKmer.get(i).r.frequency;
/*				if (!readsWhichContKmer.contains(readsOfKmer.get(i).r.nucl))
				{
					kCount1+=readsOfKmer.get(i).r.frequency;
					readsWhichContKmer.add(readsOfKmer.get(i).r.nucl);
				}*/
			}
//			int kCount = readsOfKmer.size();
			if (kCount > maxKCount)
				maxKCount = kCount;
			allKmers.add(new Kmer_general(nucl,kCount));
			for (int i = 0; i < readsOfKmer.size(); i++)
			{
				Kmer km = new Kmer(nucl,readsOfKmer.get(i).stPos);
				km.setKCount(kCount);
//				km.setKCount1(kCount1);
				readsOfKmer.get(i).r.addKmer(km);
			}
		}
			
		
		valuesDistribution = new int[maxKCount+1];
		it = kmers.entrySet().iterator();
		while (it.hasNext())
		{
			Map.Entry me = (Map.Entry)it.next();
			ArrayList<Kmer_occur> readsOfKmer = ((ArrayList<Kmer_occur>) me.getValue());
			int kCount = 0;
			for (int i = 0; i < readsOfKmer.size(); i++)
				kCount+= readsOfKmer.get(i).r.frequency;
			valuesDistribution[kCount]++;
		}
		
		valuesDistribution1 = new int[maxKCount+1];
		for (Read r : reads)
			for(Kmer km : r.kmers)
                        {
                            if (km.kCount > maxKCount)
                                DynamicOut.printStep(r.name + " " + km.StartPos);
                            valuesDistribution1[km.kCount]+=r.frequency;
                        }
		DynamicOut.finishSteps();
	}
	
	void clustering() throws IOException

        {
		 int count = 0;
                 int count1 = 0;
		 for (Read r : reads)
		 {
			 count++;
                         count1++;
                         System.out.println("Clustering: read " + count + "/" + reads.size());
			 File f = new File("kmer" + count + "data.txt");
			 FileWriter fw = new FileWriter(f);
//			 fw.write(r.kmers.size() + " " + 2 + "\n");
//			 for (Kmer km : r.kmers)
//				 fw.write(km.StartPos+ " " + km.getCount() + "\n");
			 fw.write(r.kmers.size() + " " + 1 + "\n");
			 for (Kmer km : r.kmers)
				 fw.write(km.getCount() + "\n");
			 fw.close();
			 
			 Runtime run=Runtime.getRuntime(); 
			 Process p=null;
                         if (count1 == 3000)
                         {
                            System.gc();
                            count1 = 0;
                         }
                         String exPath = "fams"+File.separator+"fams 0 0 200 " + "kmer" + count + "data ."+ File.separator;
			 p=run.exec(exPath); 
//			 p=run.exec("fams 0 0 200 " + "kmer" + count + "data .\\"); 
			 try {
				p.waitFor();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			
			 if (p.exitValue()!= 0)
			 {
				 System.out.println("Error in clustering program");
				 continue;
			 }
			 
			 File fout = new File("out_" + "kmer" + count + "data.txt");
			 
		/*	 while (!fout.exists())
			 {
				 System.out.println("Waiting when process finish");
			 }*/
			 
			 FileReader fr = new FileReader(fout);
			 BufferedReader br = new BufferedReader(fr);
			 
			 String s = br.readLine();
			 int i = 0;
			 while(s!= null)
			 {
		/*		 StringTokenizer st = new StringTokenizer(s, " \t\n");
				 st.nextToken();
				 r.kmers.get(i).setClusterParam(Double.parseDouble(st.nextToken()));*/
				 r.kmers.get(i).setClusterParam(Double.parseDouble(s));
				 i++;
				 s = br.readLine();
			 }
                         br.close();
                         fr.close();
			 
			 f.delete();
		 }
                 DynamicOut.finishSteps();
	}

        void findFreqThresholdPoisson()
        {
            double P = 1;
            // put value of lamba into Poisson class (static for simplicity)
                        int r = 1;
            while (valuesDistribution1[r] != 0)
                r++;
            double[] valuesDistributionS = new double[r-1];
            for (int i = 0; i < r-1; i++)
                valuesDistributionS[i] = valuesDistribution1[i+1];
            estimateLambdaPoisson(valuesDistributionS, 1);
            int n = 1;
            while (P > 0.01) {
                try {
                    P = 1 - pPoiss(n);
                } catch (Exception ex) {
                    Logger.getLogger(DataSet.class.getName()).log(Level.SEVERE, "", ex);
                }
                n++;
            }
            freqThr = n;

            System.out.println("Threshold =  " + n);
        }

         void findFreqThresholdExponent()
        {
            double P = 1;
            // put value of lamba into Exponent class (static for simplicity)
            int r = 1;
            while (valuesDistribution1[r] != 0)
                r++;
            double[] valuesDistributionS = new double[r-1];
            for (int i = 0; i < r-1; i++)
                valuesDistributionS[i] = valuesDistribution1[i+1];
            estimateLambdaExponent(valuesDistributionS, 1);
            int n = 1;
            while (P > 0.01) {
                try {
                    P = 1 - pExpon(n);
                } catch (Exception ex) {
                    Logger.getLogger(DataSet.class.getName()).log(Level.SEVERE, "", ex);
                }
                n++;
            }
            freqThr = n;

            System.out.println("Threshold =  " + n);
        }
        void findFreqThreshold()
        {

                if (finderrorsseglen == 0)
                {
                    int i = 1;
                    while ((valuesDistribution[i] != 0)&&(i < valuesDistribution.length))
                        i++;
                    if (i >= valuesDistribution.length)
                    {
                        System.out.println("No such number of consequtive zeros. Please specify a smaller number");
                        System.exit(1);
                    }
                    freqThr = i;
                    return;
                }

                double[] variations = new double[valuesDistribution.length+1];
                int varseglen = finderrorsseglen;
                
                if (valuesDistribution.length - varseglen <= 0)
                    {
                        System.out.println("No such number of consequtive zeros. Please specify a smaller number");
                        System.exit(1);
                    }
/*                for (int i = 0; i < valuesDistribution.length; i++)
                    System.out.println(i + " " + valuesDistribution[i]);*/
                for (int i = 1; i < valuesDistribution.length - varseglen; i++)
                {
                    double cnt = 0;
                    for (int j = i+1; j <= i + varseglen; j++)
                        cnt +=Math.abs(valuesDistribution1[j] - valuesDistribution1[j-1]);
                    variations[i] = cnt;
                }
                int i = 1;
                while ((variations[i]!= 0)&&(i < variations.length))
                    i++;
                if (i >= variations.length)
                    {
                        System.out.println("No such number of consequtive zeros. Please specify a smaller number");
                        System.exit(1);
                    }
                while ((variations[i]== 0)&&(i < variations.length))
                    i++;
                if (i >= variations.length)
                    {
                        System.out.println("No such number of consequtive zeros. Please specify a smaller number");
                        System.exit(1);
                    }
                freqThr = i+varseglen-1;
                System.out.println();
        }

	public void PrintReadsStat(String dir) throws IOException
	{
		File f = new File(dir);
		if (!f.exists())
			f.mkdir();
		int count = 1;
		int count1 = 0;
		String fol_sign = "//";
/*		for (Read r : reads)
		{
			System.out.println("Writing results: read " + count + "/" + reads.size()+"\n");	
			count++;
			count1+=r.frequency;
			String name;
			if (badReads != null)
			{
				if (badReads.contains(r))
					name = dir + fol_sign + "kmers" +"(k=" + k + ")" + "(" + r.name+")"+"bad"+"("+r.frequency+")"+ ".txt";
				else
					name = dir + fol_sign + "kmers" +"(k=" + k + ")" + "(" + r.name+")"+"good"+"("+r.frequency+")"+ ".txt";
			}
			else
				name = dir + fol_sign + "kmers" +"(k=" + k + ")" + "(" + r.name+")"+"("+r.frequency+")"+ ".txt";
			FileWriter fw = new FileWriter(name);
//			fw.write(r.kmers.size() + " " + 2 + "\n");
			int maxKcount = 0;
			for (Kmer km : r.kmers)
			{
			//	fw.write(km.nucl + " " + km.StartPos+ " " + km.getCount() + "\n");
				fw.write(km.StartPos+ " " + km.getCount() + " " + r.getNucl(km.StartPos) + " " + km.clusterParam1 + " " + r.NuclType[km.StartPos] + "\n");
				if (km.getCount() > maxKcount)
					maxKcount = km.getCount();
			}
			for (int i = r.getLength()-k+1; i < r.getLength(); i++)
				fw.write(i+ " " + "N/A" + " " + r.getNucl(i) + " " + "N/A" + " " + "N/A" + "\n");
		
			fw.close();
		}
		System.out.println("Total number of reads =" + count1);
		
		FileWriter fw = new FileWriter(dir + fol_sign + "kvalues_distribution" +"(k=" + k + ")" +"("+dir+")"+ ".txt");
		for (int i = 1; i < valuesDistribution.length; i++)
			fw.write(i + " " + valuesDistribution[i] + "\n");
		fw.close();*/



                double[] clust = new double[valuesDistribution.length+1];
                double[] clust1 = new double[valuesDistribution.length+1];
/*		File fclust = new File(file_name + "data.txt");
		FileWriter fwclust = new FileWriter(fclust);
                int z = valuesDistribution1.length-1;
                fwclust.write(z + " " + 1 + "\n");
		for (int i = 1; i < valuesDistribution1.length; i++)
                {
                         int o = valuesDistribution1[i] + i;
			 fwclust.write(o + "\n");
//                         fwclust.write(i + " " + valuesDistribution1[i] + "\n");
                }
		fwclust.close();

		Runtime run=Runtime.getRuntime();
		Process p=null;
		p=run.exec("fams//fams.exe 0 0 200 " + file_name + "data .\\");
		try {
			p.waitFor();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		if (p.exitValue()!= 0)
		{
			 System.out.println("Error in clustering program");
		}

		File fout = new File("out_" + file_name + "data.txt");


		FileReader fr = new FileReader(fout);
		BufferedReader br = new BufferedReader(fr);

		String s = br.readLine();
		int iclust = 1;
		while(s!= null)
		{
                    StringTokenizer st = new StringTokenizer(s," ");
                    while (st.hasMoreTokens())
                    {
                        clust[iclust] = Double.parseDouble(st.nextToken());
 //                       clust1[iclust] = Double.parseDouble(st.nextToken());
                    }
                    iclust++;
                    s = br.readLine();
		}

		fclust.delete();*/

                double[] variations = new double[valuesDistribution.length+1];
                int varseglen = 25;
                for (int i = 1; i < valuesDistribution.length - varseglen; i++)
                {
                    double cnt = 0;
                    for (int j = i+1; j <= i + varseglen; j++)
                        cnt +=Math.abs(valuesDistribution1[j] - valuesDistribution1[j-1]);
                    variations[i] = cnt;
                }

//                FileWriter fw1 = new FileWriter(dir + fol_sign + "kvalues_distribution1" +"(k=" + k + ")" +"("+dir+")"+ ".txt");
                FileWriter fw1 = new FileWriter(dir + fol_sign + "kvalues_distribution1.txt");
		for (int i = 1; i < valuesDistribution1.length; i++)
			fw1.write(i + " " + valuesDistribution1[i] + " " + clust[i] +" "+ clust1[i] +" " + variations[i] + "\n");
		fw1.close();


/*                FileWriter fw2 = new FileWriter(dir + fol_sign + "kvalues_distribution1_boxes" +"(k=" + k + ")" +"("+dir+")"+ ".txt");
                int count2 = 0;
                int sum = 0;
		for (int i = 1; i < valuesDistribution1.length; i++)
                {
                    if (count2 < 20)
                    {
                        count2++;
                        sum+=valuesDistribution1[i];
                    }
                    else
                    {
			fw2.write(i + " " + sum + "\n");
                        sum = 0;
                        count2 = 0;
                    }
                }
		fw2.close();*/
	}
	void clear()
	{
		int count = 0;
		 for (Read r : reads)
		 {
			 count++;
//			 System.out.println("Deleting files: read " + count + "/" + reads.size()+"\n");	
			 File f = new File("modes_kmer" + count + "data.txt");			
			 f.deleteOnExit();
			 f = new File("pilot_200_kmer" + count + "data.txt");
			 f.deleteOnExit();
			 f = new File("out_kmer" + count + "data.txt");
                         f.delete();
			 int i;
			 f.deleteOnExit();
			 
		 }
	}
	
	void deleteStrangeNucl()
	{
		ArrayList<Read> toDel = new ArrayList<Read>();
		for (Read r : reads)
		{
			if (r.nucl.contains("N"))
				toDel.add(r);
		}
		reads.removeAll(toDel);
			
	}
	void findErrors()
	{
		badReads = new ArrayList<Read>();
		goodReads = new ArrayList<Read>();
		for (Read r : reads)
		{
			boolean bad = false;
			for (Kmer k : r.kmers)
				if (k.kCount <= freqThr)
				{
					badReads.add(r);
					bad = true;
					break;
				}
			if (!bad)
				goodReads.add(r);
		}
		int bad = 0;
		int all = 0;
		for (Read r : badReads)
			bad+=r.frequency;
		for (Read r : reads)
			all+=r.frequency;
		System.out.println("bad = " + bad + "/ all = " + all);
		
/*		HashSet hs = new HashSet();
		for (Read r : goodReads)
			for (Kmer k: r.kmers)
				hs.add(k);
		
		System.out.println(hs.size());*/
		
		for (Read r : badReads)
		{
//			r.errorRegions = new ArrayList<ErrorRegion>();
			int a[] = new int[r.getLength()-k+1];
			double b[] = new double[r.getLength()-k+1];
			for (Kmer k : r.kmers)
			{
				a[k.getStartPos()] = k.getKCount();
				b[k.getStartPos()] = k.clusterParam1;
			}
			int i = 0;
			int l = 0;
			String s = "";
			int start = 0;
			double clustprm = 0;
			while (i <= r.getLength()-k)
			{
				if (a[i] <= freqThr)
				{
					if (l==0)
					{
						start = i;
						clustprm = b[i];
					}
					l++;
					s+=r.getNucl(i);
					r.NuclType[i] = 1;
				}
				else
				{
					if (l > 0)
					{
						if ((Math.abs(b[i] - clustprm) < 0.1) && (a[i] <= freqThr*additionalfreqThrMult)&&(clustprm >= 0))
						{
							l++;
							s+=r.getNucl(i);
							r.NuclType[i] = 1;
						}
						else
						{
							int upd_start = start-1;
							while (upd_start >=0)
							{
								if ((Math.abs(b[upd_start] - clustprm) < 0.1) && (a[upd_start] <= freqThr*additionalfreqThrMult)&&(clustprm >= 0))
								{
									l++;
									s = r.getNucl(upd_start) + s;
									r.NuclType[upd_start] = 1;
									start = upd_start;
									upd_start--;
								}
								else
									break;
							}
							r.addErrorRegion(new ErrorRegion(s,start,i-1,l));
							l = 0;
							s = "";
						}
					}
					else
					{
						l = 0;
						s = "";
					}
				}
				i++;
			}
			if (l > 0)
			{
				int upd_start = start-1;
				while (upd_start >=0)
				{
					if ((Math.abs(b[upd_start] - clustprm) < 0.1) && (a[upd_start] <= freqThr*additionalfreqThrMult)&&(clustprm >= 0))
					{
						l++;
						s = r.getNucl(upd_start) + s;
						r.NuclType[upd_start] = 1;
						start = upd_start;
						upd_start--;
					}
					else
						break;
				}
				r.addErrorRegion(new ErrorRegion(s,start,r.getLength()-k,l));
			}
			
			i = 0;
/*			while (i < r.errorRegions.size()-1)
			{
				if (r.errorRegions.get(i).end >= r.errorRegions.get(i+1).begin)
				{
// slow, can be faster					
					String union = r.nucl.substring(r.errorRegions.get(i).begin, r.errorRegions.get(i+1).end+1);
					int s_new = r.errorRegions.get(i).begin;
					int e_new = r.errorRegions.get(i+1).end;
					ErrorRegion r_new = new ErrorRegion(union,s_new,e_new,e_new-k+1-s_new);
					r.errorRegions.remove(i+1);
					r.errorRegions.remove(i);
					r.errorRegions.add(i, r_new);
				}
				else
					i++;
			}*/
		}
	}
	
	void printErrorsStat(String dir) throws IOException
	{
		File f = new File(dir);
		if (!f.exists())
			f.mkdir();
		String fol_sign = "//";
		String name = dir + fol_sign + "error_stat(k=" + k +").txt";
		FileWriter fw = new FileWriter(name); 
		
// all errors
		for (Read r : badReads)
			for (ErrorRegion er : r.errorRegions)
			{
				String type = new String();
				if (er.begin == 0)
					type = "begin";
				if((er.end + k  == r.getLength())&&(er.begin != 0))
					type = "end";
				if((er.end + k  != r.getLength())&&(er.begin != 0))
				{
					if (r.getNucl(er.end) == r.getNucl(er.end+1))
						type = "ins";
					else
						type = "delOrRepl";
				}
				
				String s = new String();
				s = er.length + " " + r.name  + " " + er.begin + " " + type + " " + (double) er.length/r.getLength()+ " " + r.getLength()+ "\n";
				fw.write(s);
			}
		
		fw.close();
		
		// how many errors of different lengths
		fw = new FileWriter(dir+fol_sign +"error_stat_lendistrib_(k=" + k +").txt");
		int ml = 0;
		for (Read r : badReads)
			for (ErrorRegion er : r.errorRegions)
				if (er.length > ml)
					ml = er.length;
		int[] a = new int[ml+1];
		for (Read r : badReads)
		{
			for (ErrorRegion er : r.errorRegions)
				a[er.length]++;
		}
		for (int i = 1; i <= ml; i++)
			fw.write(i + " " + a[i]+"\n");
		fw.close();
		
		FileWriter fw1 = new FileWriter(dir+fol_sign +"kerrors_number_repl(k=" + k +").txt");
		FileWriter fw2 = new FileWriter(dir+fol_sign +"(k-2)errors_stat(k=" + k +").txt");
		FileWriter fw3 = new FileWriter(dir+fol_sign +"(k-3)errors_stat(k=" + k +").txt");
		FileWriter fw4 = new FileWriter(dir+fol_sign +"(k-4)errors_stat(k=" + k +").txt");
		FileWriter fwm1 = new FileWriter(dir+fol_sign +"(k-1)errors_stat(k=" + k +").txt");
		fw = new FileWriter(dir+fol_sign +"kerrors_stat(k=" + k +").txt");
		for (Read r : badReads)
			for (ErrorRegion er : r.errorRegions)
			{
				if ((er.length == k)&&(er.begin!= 0)&&(er.end + k != r.getLength()))
				{
						String s = r.nucl.substring(er.begin, er.end);
						char c = r.nucl.charAt(er.end);
						char repl = c;
						int maxKcount = 0;
                                                int count = 0;
						for (Kmer_general km : allKmers)
							if ((km.kCount > freqThr)&&(km.nucl.startsWith(s)))
								if (km.nucl.charAt(k-1)!= c)
                                                                {
                                                                        count++;
									if (km.kCount > maxKcount)
									{
										maxKcount = km.kCount;
                                                                                repl = km.nucl.charAt(k-1);
									}
                                                                }
						fw.write(r.name + " " + er.begin + " " + maxKcount + " " + count + " ");
                                                if (repl == r.getNucl(er.end+1))
                                                    fw.write("ins" + "\n");
                                                else
                                                    fw.write("repl" + "\n");
					
				}
				if (er.length == k-1)
				{
//					System.out.println("(k-1)-errors considering: " + r.name);
					fwm1.write(r.name  + " " + er.begin + " ");
					if (r.getNucl(er.end) == r.getNucl(er.end+1))
						fwm1.write(1 + "\n");
					else
						fwm1.write(0 + "\n");
				}
				if (er.length == k-2)
				{
//					System.out.println("(k-2)-errors considering: " + r.name);
					fw2.write(r.name  + " " + er.begin + " ");
					int l = 0;
					int i = er.end+1;
					while((r.getNucl(i) == r.getNucl(er.end))&&(i < r.getLength()))
					{
						l++;
						i++;
					}
					fw2.write(l + " ");
					String s = r.getNucl(er.begin, er.end+1);
					s = s + r.getNucl(er.end+1);
					int bestKcount = 0;
					for (Kmer_general km : allKmers)
						if ((km.kCount > freqThr)&&(km.nucl.contentEquals(s)))
						{
							if (km.kCount > bestKcount)
								bestKcount = km.kCount;
						}
					fw2.write(bestKcount + "\n");
					
				}
				if (er.length == k-3)
				{
//					System.out.println("(k-3)-errors considering: " + r.name);
					fw3.write(r.name  + " " + er.begin + " ");
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
						while((r.getNucl(j) == r.getNucl(er.end+1))&&(j < r.getLength()))
						{
							q++;
							j++;
						}
					fw3.write(l + " " + q + " ");
					int bestKcount = 0;
					if (q == 1)
					{
						String s = r.getNucl(er.begin, er.end+1+q);
						s = s + r.getNucl(er.end+1);
						
						for (Kmer_general km : allKmers)
							if ((km.kCount > freqThr)&&(km.nucl.contentEquals(s)))
							{
								if (km.kCount > bestKcount)
									bestKcount = km.kCount;
							}
					}
					fw3.write(bestKcount + "\n");
					
				}
				if (er.length == k-4)
				{
//					System.out.println("(k-4)-errors considering: " + r.name);
					fw4.write(r.name  + " " + er.begin + " ");
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
						while((r.getNucl(j) == r.getNucl(er.end+1))&&(j < r.getLength()))
						{
							q++;
							j++;
						}
					fw4.write(l + " " + q + " ");
					int bestKcount = 0;
					if (q == 2)
					{
						String s = r.getNucl(er.begin, er.end+1+q);
						s = s + r.getNucl(er.end+1);
						
						for (Kmer_general km : allKmers)
							if ((km.kCount > freqThr)&&(km.nucl.contentEquals(s)))
							{
								if (km.kCount > bestKcount)
									bestKcount = km.kCount;
							}
					}
					fw4.write(bestKcount + "\n");
					
				}
			}
		
		System.out.println(allKmers.size());
		int bad = 0;
		int good = 0;
		for (Read r : badReads)
			bad += r.frequency;
		for (Read r : goodReads)
			good += r.frequency;
		System.out.println(bad + "/" + good);
		fw.close();
		fw1.close();
		fw2.close();
		fw3.close();
		fw4.close();
		fwm1.close();
		
		FileWriter fws = new FileWriter(dir+fol_sign +"strange_nucl(k=" + k +").txt");
		for (Read r : reads)
		{
			for (int i = 0; i < r.getLength(); i++)
				if ((r.getNucl(i)!='A')&&(r.getNucl(i)!='T')&&(r.getNucl(i)!='G')&&(r.getNucl(i)!='C'))
					fws.write(r.name + " " + i + "\n");
		}
		fws.close();

/*		int i = 0;
		for (Read r : reads)
		{
			if(r.nucl.contains("GACCAGTTCATCATCATATCCCTAN"))
				i++;
		}
		System.out.println("AAAAA=" + i);*/
		FileWriter fw_er_len_ratio = new FileWriter(dir+fol_sign +"errors_reg_len_ratio(k=" + k +").txt");
		for (Read r : badReads)
		{
			double l = 0;
			int lr = r.getLength();
			for (ErrorRegion er : r.errorRegions)
				l+=er.length;
			double d =  100.0*(double)(l/lr);
			fw_er_len_ratio.write(r.name + " " + l + " " + r.getLength() + " " + d + "\n");
		}
		fw_er_len_ratio.close();
		
	}
	void cutTails()
	{
		for (Read r : badReads)
		{
			System.out.println("Deleting tails: " + r.name +"\n" );
			ArrayList<ErrorRegion> corrected = new ArrayList<ErrorRegion>();
			for (ErrorRegion er : r.errorRegions)
			{
				if (er.begin == 0)
					for (int i = er.begin; i <= er.end; i++)
						r.corrections[i] = new Correction('D');
				
				if ((er.end + k  == r.getLength())&&(er.begin != 0))
					for (int i = er.begin+k-1; i < r.getLength(); i++)
						r.corrections[i] = new Correction('D');
			/*	if (er.length > r.getLength()*(maxErrorlenPerc/100))
					reads.remove(r);*/
			}
		}
	}
	void correctErrors_k()
	{
		
		for (Read r : badReads)
			for (ErrorRegion er : r.errorRegions)
				if ((er.length == k)&&(er.begin!= 0)&&(er.end + k != r.getLength()))
				{
					System.out.println("k-errors correction: " + r.name);
					String s = r.nucl.substring(er.begin, er.end);
					char c = r.nucl.charAt(er.end);
					char repl = c;
					int maxKcount = 0;
					for (Kmer_general km : allKmers)
						if ((km.kCount > freqThr)&&(km.nucl.startsWith(s)))
							if (km.nucl.charAt(k-1)!= c)
								if (km.kCount > maxKcount)
								{
									maxKcount = km.kCount;
									repl = km.nucl.charAt(k-1);
								}
					r.corrections[er.end] = new Correction('R',""+repl);
				}
		
		
	}
	void correctErrors_kminus1()
	{
		
		for (Read r : badReads)
			for (ErrorRegion er : r.errorRegions)
				if ((er.length == k-1)&&(er.begin!= 0)&&(er.end + k < r.getLength()-1))
				{
					System.out.println("(k-1)-errors correction: " + r.name);
					if (r.getNucl(er.end) == r.getNucl(er.end+1))
						r.corrections[er.end] = new Correction('D');
				}
		
		
	}
	void correctErrors_kminus2()
	{
		
		for (Read r : badReads)
			for (ErrorRegion er : r.errorRegions)
				if (er.length == k-2)
				{
					int l = 0;
					int i = er.end+1;
					while((r.getNucl(i) == r.getNucl(er.end))&&(i < r.getLength()))
					{
						l++;
						i++;
					}
					String s = r.getNucl(er.begin, er.end+1);
					s = s + r.getNucl(er.end+1);
					int bestKcount = 0;
					for (Kmer_general km : allKmers)
						if ((km.kCount > freqThr)&&(km.nucl.contentEquals(s)))
						{
							if (km.kCount > bestKcount)
								bestKcount = km.kCount;
						}
					if (l == 2)
						r.corrections[er.end] = new Correction('D');
					if (l==0)
						if (bestKcount > 0)
							r.corrections[er.end+1] = new Correction('I',r.getNucl(er.end+1));
				}
		
		
	}
	void correctErrors_kminus3()
	{
		
		for (Read r : badReads)
			for (ErrorRegion er : r.errorRegions)
				if (er.length == k-3)
				{
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
						while((r.getNucl(j) == r.getNucl(er.end+1))&&(j < r.getLength()))
						{
							q++;
							j++;
						}
					int bestKcount = 0;
					if (q == 1)
					{
						String s = r.getNucl(er.begin, er.end+1+q);
						s = s + r.getNucl(er.end+1);
						
						for (Kmer_general km : allKmers)
							if ((km.kCount > freqThr)&&(km.nucl.contentEquals(s)))
							{
								if (km.kCount > bestKcount)
									bestKcount = km.kCount;
							}
					}
					if (l == 3)
						r.corrections[er.end] = new Correction('D');
					if (q==1)
						if (bestKcount > 0)
							r.corrections[er.end+1+q] = new Correction('I',r.getNucl(er.end+1));
				}
		
		
	}
	void correctErrors_kminus4()
	{		
		for (Read r : badReads)
			for (ErrorRegion er : r.errorRegions)
				if (er.length == k-4)
				{
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
						while((r.getNucl(j) == r.getNucl(er.end+1))&&(j < r.getLength()))
						{
							q++;
							j++;
						}
					int bestKcount = 0;
					if (q == 2)
					{
						String s = r.getNucl(er.begin, er.end+1+q);
						s = s + r.getNucl(er.end+1);
						
						for (Kmer_general km : allKmers)
							if ((km.kCount > freqThr)&&(km.nucl.contentEquals(s)))
							{
								if (km.kCount > bestKcount)
									bestKcount = km.kCount;
							}
					}
					if (l == 4)
						r.corrections[er.end] = new Correction('D');
					if (q==2)
						if (bestKcount > 0)
							r.corrections[er.end+1+q] = new Correction('I',r.getNucl(er.end+1));
				}	
	}
	void PrintCorrectedReads(String outfile) throws IOException

	{
		FileWriter fw = new FileWriter(outfile);
		for(Read r : reads)
		{
			String s = "";
			for (int i = 0; i < r.getLength(); i++)
				if (r.corrections[i] != null)
				{
					if (r.corrections[i].type == 'R')
						s = s + r.corrections[i].replacement;
					if (r.corrections[i].type == 'I')
						s = s + r.nucl.charAt(i) + r.corrections[i].replacement;
				}
				else
					s = s + r.nucl.charAt(i);
			if (s.length() > lengthThr)
				for (int i = 0; i < r.frequency; i++)
					fw.write(">" + r.name + "_" + i + "\n" + s + "\n");
		}
		fw.close();
	}	
	public void PrintHaplotypes(String outfile) throws IOException
	{
		FileWriter fw = new FileWriter(outfile);
		HaplComp hc = new HaplComp();
		Collections.sort(haplotypes, hc);
		
		double count = 0;
		for (Haplotype h : haplotypes)
			count += h.concentration;
		for (Haplotype h : haplotypes)
		{
			int decimalPlace = 4;
			BigDecimal bd = new BigDecimal(100*h.concentration);
			bd = bd.setScale(decimalPlace,BigDecimal.ROUND_UP);
			double x = bd.doubleValue();
			bd = new BigDecimal(100*h.concentration/count);
			bd = bd.setScale(decimalPlace,BigDecimal.ROUND_UP);
			double y = bd.doubleValue();
			fw.write(">" + h.name + "_" + h.frequency + "_" + x + "_"+ "|" + y + "|" + "\n" + h.nucl + "\n");
		}
		fw.close();
	}
        public void PrintHaplotypesWithNameTag(String outfile, String tag) throws IOException
	{
		FileWriter fw = new FileWriter(outfile);
		HaplComp hc = new HaplComp();
		Collections.sort(haplotypes, hc);
		
		double count = 0;
		for (Haplotype h : haplotypes)
			count += h.concentration;
		for (Haplotype h : haplotypes)
		{
			int decimalPlace = 4;
			BigDecimal bd = new BigDecimal(100*h.concentration);
			bd = bd.setScale(decimalPlace,BigDecimal.ROUND_UP);
			double x = bd.doubleValue();
			bd = new BigDecimal(100*h.concentration/count);
			bd = bd.setScale(decimalPlace,BigDecimal.ROUND_UP);
			double y = bd.doubleValue();
			fw.write(">"+ h.name + "_" + tag + "_" + h.frequency + "_" + x + "_"+ "|" + y + "|" + "\n" + h.nucl + "\n");
		}
		fw.close();
	}
	public void PrintUniqueReads(String outfile) throws IOException
	{
		FileWriter fw = new FileWriter(outfile);
		for (Read r : reads)
		{
			fw.write(">" + r.name + "_" + r.frequency + "\n" + r.nucl + "\n");		
		}
		fw.close();
	}
	public void findHaplotypes() throws IOException
	{
		haplotypes = new ArrayList<Haplotype>();
		StrComp sc = new StrComp();
		Collections.sort(reads, sc);
		System.out.println("Finding haplotypes");
                int goodhapl = 0;
                int badhapl = 0;
                int lenratio = 90; //90
		for (Read r : reads)
		{
			boolean toAdd = true;
			int nhap = 0;
			for (Haplotype h : haplotypes)
				if (h.containsRead(r))
				{
					toAdd = false;
					nhap++;
				}
			if (toAdd)
                        {
				haplotypes.add(new Haplotype(r));
                       /*         if (badReads.contains(r))
                                    badhapl++;
                                else
                                    goodhapl++;*/
                        }
			else
				for (Haplotype h : haplotypes)
					if (h.containsRead(r))
						h.concentration += r.frequency/nhap;
		}
		for (Haplotype h : haplotypes)
			h.concentration = h.concentration / getNreads();


                ArrayList<Haplotype> shorthapl = new ArrayList<Haplotype>();
                double allconc = 0;
                for (Haplotype h : haplotypes)
                    allconc+=h.concentration;
                double avlen = 0;
                for (Haplotype h : haplotypes)
                    avlen += h.getLength()*((double) h.concentration/allconc);
                for (Haplotype h : haplotypes)
                    if (100 * ((double) h.getLength()/avlen) < lenratio)
                        shorthapl.add(h);
                haplotypes.removeAll(shorthapl);
                for (Haplotype h : haplotypes)
                    allconc+=h.concentration;
                for (Haplotype h : haplotypes)
                    h.concentration /= allconc;


                System.out.println("Goodhapl=" + goodhapl);
                System.out.println("Badhapl=" + badhapl);
	}
	void compareHaplWithClones(String clonefile) throws IOException
	{
		DataSet cds = new DataSet(clonefile,'c');
		for (Haplotype h : haplotypes)
		{
			String clns = "";
			for (Read r : cds.reads)
				if (h.nucl.contains(r.nucl) || r.nucl.contains(h.nucl))
					clns = clns + r.name;
			if (clns != "")
				h.clone = clns;
		}
	}
	void delLargeErrors()
	{
		for (Read r : badReads)
		{
			double l = 0;
			int lr = r.getLength();
			for (ErrorRegion er : r.errorRegions)
				l+=er.length;
			double d =  100.0*(double)(l/lr);
			if (d > r.getLength()*(maxAllErrorsPerc/100))
				reads.remove(r);
		}
	}
	int getNerrors()
	{
		int count = 0;
		for (Read r : reads)
			count+= r.errorRegions.size();
		return count;
	}
	int maxReadLen()
	{
		int m = 0;
		for (Read r : reads)
			if (r.getLength() > m)
				m = r.getLength();
		return m;
	}
	int minReadLen()
	{
		int m = Integer.MAX_VALUE;
		for (Read r : reads)
			if (r.getLength() < m)
				m = r.getLength();
		return m;
	}
	int getNreads()
	{
		int n = 0;
		for (Read r : reads)
			n+= r.frequency;
		return n;
	}
	int getNuniquereads()
	{
		return reads.size();
	}
	double AvReadLen()
	{
		double m = 0;
		for (Read r : reads)
			m += r.getLength();
		return m/reads.size();
	}
	int getNUniqueErrorReads()
	{
		return badReads.size();
	}
        void PrintReadsDistrib(String dir) throws IOException
	{
		String fol_sign = "//";

		FileWriter fw = new FileWriter(dir + fol_sign + "kvalues_distribution" +"(k=" + k + ")" +"("+dir+")"+ ".txt");
		for (int i = 1; i < valuesDistribution.length; i++)
			fw.write(i + " " + valuesDistribution[i] + "\n");
		fw.close();

                double[] clust = new double[valuesDistribution.length+1];
                double[] clust1 = new double[valuesDistribution.length+1];
/*		File fclust = new File(file_name + "data.txt");
		FileWriter fwclust = new FileWriter(fclust);
                int z = valuesDistribution1.length-1;
                fwclust.write(z + " " + 1 + "\n");
		for (int i = 1; i < valuesDistribution1.length; i++)
                {
                         int o = valuesDistribution1[i] + i;
			 fwclust.write(o + "\n");
//                         fwclust.write(i + " " + valuesDistribution1[i] + "\n");
                }
		fwclust.close();

		Runtime run=Runtime.getRuntime();
		Process p=null;
		p=run.exec("fams//fams.exe 0 0 200 " + file_name + "data .\\");
		try {
			p.waitFor();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		if (p.exitValue()!= 0)
		{
			 System.out.println("Error in clustering program");
		}

		File fout = new File("out_" + file_name + "data.txt");


		FileReader fr = new FileReader(fout);
		BufferedReader br = new BufferedReader(fr);

		String s = br.readLine();
		int iclust = 1;
		while(s!= null)
		{
                    StringTokenizer st = new StringTokenizer(s," ");
                    while (st.hasMoreTokens())
                    {
                        clust[iclust] = Double.parseDouble(st.nextToken());
 //                       clust1[iclust] = Double.parseDouble(st.nextToken());
                    }
                    iclust++;
                    s = br.readLine();
		}

		fclust.delete();*/

                double[] variations = new double[valuesDistribution.length+1];
                int varseglen = 25;
                for (int i = 1; i < valuesDistribution.length - varseglen; i++)
                {
                    double cnt = 0;
                    for (int j = i+1; j <= i + varseglen; j++)
                        cnt +=Math.abs(valuesDistribution1[j] - valuesDistribution1[j-1]);
                    variations[i] = cnt;
                }

                FileWriter fw1 = new FileWriter(dir + fol_sign + "kvalues_distribution1" +"(k=" + k + ")" +"("+dir+")"+ ".txt");
		for (int i = 1; i < valuesDistribution1.length; i++)
			fw1.write(i + " " + valuesDistribution1[i] + " " + clust[i] +" "+ clust1[i] +" " + variations[i] + "\n");
		fw1.close();


/*                FileWriter fw2 = new FileWriter(dir + fol_sign + "kvalues_distribution1_boxes" +"(k=" + k + ")" +"("+dir+")"+ ".txt");
                int count2 = 0;
                int sum = 0;
		for (int i = 1; i < valuesDistribution1.length; i++)
                {
                    if (count2 < 20)
                    {
                        count2++;
                        sum+=valuesDistribution1[i];
                    }
                    else
                    {
			fw2.write(i + " " + sum + "\n");
                        sum = 0;
                        count2 = 0;
                    }
                }
		fw2.close();*/
	}
        void PrintReadsFreqDistrib(String addr) throws IOException
        {
            FileWriter fw = new FileWriter(addr);
            for (Read r:reads)
                fw.write(r.frequency + "\n");
            fw.close();
        }
        void calcDirection(String ref_file, double gapop, double gapext) throws IOException
        {
            DataSet ref = new DataSet(ref_file);
            FileWriter fw = new FileWriter(this.file_name +"_refStat.txt");
            for (Read r : reads)
            {
                fw.write(r.name + " ");
                for (int i = 0; i < 2; i++)
                {
                    String s_ref = ref.reads.get(0).nucl;
                    String q = "Direct";
                    Runtime run=Runtime.getRuntime();
                    Process p=null;
                    FileWriter fw_alin = new FileWriter("allign_input.fas");
                    fw_alin.write(r.name + "\n");
                    fw_alin.write(r.nucl + "\n");
                    fw_alin.write(">reference" + "\n");
                    fw_alin.write(s_ref + "\n");
                    fw_alin.close();
                    String param = "ClustalW2//clustalw2 -INFILE=" +"allign_input.fas" + " -OUTFILE=" + "allign_output.fas";
                    param+= " -OUTPUT=FASTA -DNAMATRIX=IUB -GAPOPEN=" + gapop +" -GAPEXT=" + gapext + " -TYPE=DNA -PWDNAMATRIX=IUB -PWGAPOPEN=" + gapop+ " -PWGAPEXT=" + gapext;

                    p = run.exec(param);
                    BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
                    String s1;
                    while ((s1 = stdInput.readLine()) != null) {
    //                  System.out.println(s1);
                    }
                    if (p.exitValue()!= 0)
                    {
                        System.out.println("Error in allgnment program");
                    }
                    DataSet allignment = new DataSet("allign_output.fas",'c');
                    String al_read = "";
                    String al_ref = "";
                    for (Read r1 : allignment.reads)
                        if (!r1.name.contains("reference"))
                            al_read = r1.getNucl();
                        else
                            al_ref = r1.getNucl();
                    int match = 0;
                    for (int j = 0; j < al_read.length(); j++)
                        if (al_read.charAt(j) == al_ref.charAt(j))
                            match++;
                    double score = (double) match/r.getLength();
                    fw.write(q + " " + score + " ");
                    String rev_ref = "";
                    for (int j = s_ref.length() - 1; j >=0; j --)
                    {
                        if (s_ref.charAt(j) == 'A')
                            rev_ref += 'T';
                        if (s_ref.charAt(j) == 'T')
                            rev_ref += 'A';
                        if (s_ref.charAt(j) == 'G')
                            rev_ref += 'C';
                        if (s_ref.charAt(j) == 'C')
                            rev_ref += 'G';
                    }
                    s_ref = rev_ref;
                    q = "Reverse";
                    File f = new File("allign_input.fas");
                    f.delete();
                    f = new File("allign_output.fas");
                    f.delete();
                    f = new File("allign_input.dnd");
                    f.delete();
                }
                fw.write("\n");
            }
            fw.close();
        }
        void convertZoyaData(String addr) throws FileNotFoundException, IOException
        {
                FileReader fr = new FileReader(addr);
		BufferedReader br = new BufferedReader(fr);
                FileWriter fw = new FileWriter(addr + "_converted.fas");
	
		String s = br.readLine();
                while(s.equalsIgnoreCase(""))
                    s = br.readLine();

		String nucl = new String();
		String name = new String();
		name = "";
		int i = 1;
		int j = s.length();
		while ((i < s.length())&&(s.charAt(i)!=' '))
		{
			name+=s.charAt(i);
			i++;
		}
		s = br.readLine();
		try
		{
		while (s!= null)
		{
//			System.out.println("Reading " + name);

			if (s.length() == 0)
			{
				s = br.readLine();
				continue;
			}
			if(s.charAt(0) == '>')
                        {
                                        StringTokenizer st = new StringTokenizer(name,"_");
                                        st.nextToken();
                                        st.nextToken();
                                        st.nextToken();
                                        st.nextToken();
                                        st.nextToken();
                                        int mult = Integer.parseInt(st.nextToken());
                                        for (int u = 0; u < mult; u++)
                                            fw.write(">" + name + "_" + u + "\n" + nucl + "\n");
					nucl = "";
					name = "";
					i = 1;
					while ((i < s.length())&&(s.charAt(i)!=' '))
					{
						name+=s.charAt(i);
						i++;
					}
			}
			else
				nucl+=s.toUpperCase();
			s = br.readLine();
		}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	//	reads.add(new Read(nucl, name));
                StringTokenizer st = new StringTokenizer(name,"_");
                st.nextToken();
                st.nextToken();
                st.nextToken();
                st.nextToken();
                st.nextToken();
                int mult = Integer.parseInt(st.nextToken());
                for (int u = 0; u < mult; u++)
                    fw.write(">" + name + "_" + u + "\n" + nucl + "\n");
		fr.close();
                fw.close();
        }
        void calculateRealErrorsStat(String ref, int homoplen, double gapop, double gapext) throws IOException
        {
            Read clone = (new DataSet(ref)).reads.get(0);
            clone.name = "Clone";
            FileWriter fw = new FileWriter(this.file_name + "_" + ref + "_realerrors.txt");
            fw.write("read freq len averrorsdist repl ins del ");
            for (int i = 1; i <= homoplen; i++)
                fw.write("homop" + i + " ");
            fw.write("allerr nreplat nreplac nreplag nrepltc nrepltg nreplgc" + "\n");
            FileWriter fw_un = new FileWriter(this.file_name + "_" + ref + "_realerrorsun.txt");
            fw_un.write("read freq len averrorsdist repl ins del ");
            for (int i = 1; i <= homoplen; i++)
                fw_un.write("homop" + i + " ");
            fw_un.write("allerr nreplat nreplac nreplag nrepltc nrepltg nreplgc" + "\n");
            
            
            for (Read r : reads)
            {
//                    DynamicOut.printStep(r.name + "//" + reads.size());
                    System.out.println(r.name + "//" + reads.size());
                    Runtime run=Runtime.getRuntime();
                    Process p=null;
                    FileWriter fw_alin = new FileWriter("allign_input.fas");
                    String ncl = "";
                    for (int i = r.nucl.length()-1; i>=0; i--)
                    {
                        if (r.nucl.charAt(i) == 'A')
                               ncl+= "T";
                        if (r.nucl.charAt(i) == 'T')
                               ncl+="A";
                        if (r.nucl.charAt(i) == 'G')
                               ncl+="C";
                        if (r.nucl.charAt(i) == 'C')
                               ncl+="G";
                    }
                    fw_alin.write(">" + r.name + "\n" + ncl + "\n");
                    fw_alin.write(">" + clone.name + "\n" + clone.nucl + "\n");
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
/*                    File f = new File("allign_input.fas");
                    f.delete();
                    f = new File("allign_output.fas");
                    f.delete();
                    f = new File("allign_input.dnd");
                    f.delete();*/
                    String read_align = "";
                    String clone_align = "";
                    for (Read r1: alignment.reads)
                    {
                        if (r1.name.equalsIgnoreCase(r.name))
                            read_align = r1.nucl;
                        if (r1.name.equalsIgnoreCase(clone.name))
                            clone_align = r1.nucl;
                    }
                    
                    int l = read_align.length();
                    int pos_r = -1;
                    int pos_c = -1;

                    int startcompdiff = 0;
                    int endcompdiff = l-1;
                    String sz = "";
                    if (read_align.charAt(0) == '-')
                        sz = read_align;
                    else
                        sz = clone_align;
                    while (sz.charAt(startcompdiff) == '-')
                        startcompdiff++;
                    if (read_align.charAt(l-1) == '-')
                        sz = read_align;
                    else
                        sz = clone_align;
                    while (sz.charAt(endcompdiff) == '-')
                        endcompdiff--;

                    int errors[] = new int[homoplen+2];
                    int nrepl = 0;
                    int nins = 0;
                    int ndel = 0;
                    int nreplat = 0;
                    int nreplac = 0;
                    int nreplag = 0;
                    int nrepltc = 0;
                    int nrepltg = 0;
                    int nreplgc = 0;
                    ArrayList<Integer> erpos = new ArrayList<Integer>();
                    for (int i = startcompdiff; i <= endcompdiff; i++)
                    {
                        if (read_align.charAt(i) == clone_align.charAt(i))
                            continue;
                        if ((read_align.charAt(i)!= '-')&&(clone_align.charAt(i)!='-'))
                        {
                            nrepl++;
                            int readpos = 0;
                            for (int j = 0; j <= i; j++)
                                if (read_align.charAt(j)!='-')
                                    readpos++;
                            erpos.add(readpos);
                            if ((read_align.charAt(i) == 'A')&&(clone_align.charAt(i) == 'T'))
                                nreplat++;
                            if ((read_align.charAt(i) == 'T')&&(clone_align.charAt(i) == 'A'))
                                nreplat++;
                            if ((read_align.charAt(i) == 'A')&&(clone_align.charAt(i) == 'C'))
                                nreplac++;
                            if ((read_align.charAt(i) == 'C')&&(clone_align.charAt(i) == 'A'))
                                nreplac++;
                            if ((read_align.charAt(i) == 'A')&&(clone_align.charAt(i) == 'G'))
                                nreplag++;
                            if ((read_align.charAt(i) == 'G')&&(clone_align.charAt(i) == 'A'))
                                nreplag++;
                            if ((read_align.charAt(i) == 'T')&&(clone_align.charAt(i) == 'C'))
                                nrepltc++;
                            if ((read_align.charAt(i) == 'C')&&(clone_align.charAt(i) == 'T'))
                                nrepltc++;
                            if ((read_align.charAt(i) == 'T')&&(clone_align.charAt(i) == 'G'))
                                nrepltg++;
                            if ((read_align.charAt(i) == 'G')&&(clone_align.charAt(i) == 'T'))
                                nrepltg++;
                            if ((read_align.charAt(i) == 'G')&&(clone_align.charAt(i) == 'C'))
                                nreplgc++;
                            if ((read_align.charAt(i) == 'C')&&(clone_align.charAt(i) == 'G'))
                                nreplgc++;
                            continue;
                        }
                        char z;
                        if (read_align.charAt(i)!= '-')
                            z = read_align.charAt(i);
                        else
                            z = clone_align.charAt(i);
                        int bhomr = i;
                        int ehomr = i;
                        int bhomc = i;
                        int ehomc = i;
                        while (((read_align.charAt(ehomr)==z)||(read_align.charAt(ehomr) == '-'))&&(ehomr <= endcompdiff))
                            ehomr++;
                        while ((read_align.charAt(bhomr)==z)||(read_align.charAt(bhomr) == '-')&&(bhomr >= startcompdiff))
                            bhomr--;
                        while (((clone_align.charAt(ehomc)==z)||(clone_align.charAt(ehomc) == '-'))&&(ehomc <= endcompdiff))
                            ehomc++;
                        while ((clone_align.charAt(bhomc)==z)||(clone_align.charAt(bhomc) == '-')&&(bhomc >= startcompdiff))
                            bhomc--;
                        int b = Math.max(bhomr+1, bhomc+1);
                        int e = Math.min(ehomr-1, ehomc-1);
                        if (e == b)
                        {
                            if (read_align.charAt(i) == '-')
                                ndel++;
                            else
                                nins++;
                        }
                        else
                        {
                            int len = 0;
                            for (int j = b; j <= e; j++)
                                if (clone_align.charAt(j) == z)
                                    len++;
                            errors[len]++;
                        }
                        int q = b;
                        for (int j = b; j <=e; j++)
                            if (read_align.charAt(q)==z)
                            {
                                q = j;
                                break;
                            }
                        int readpos = 0;
                        for (int j = 0; j <= q; j++)
                            if (read_align.charAt(j)!='-')
                                readpos++;
                        erpos.add(readpos);
//                        int len = Math.min(e-b+1, homoplen);
                        i = e;
                    }
                    int erdist = 0;
                    for (int i = 1; i < erpos.size(); i++)
                        erdist += (erpos.get(i)-erpos.get(i-1)-1);
                    double averrorsdist = -1;
                    if (erpos.size() >= 2)
                        averrorsdist = (double) erdist/(erpos.size()-1);
                    int allerr = nrepl + nins + ndel;
                    for (int i = 1; i <= homoplen; i++)
                        allerr += errors[i];
                    fw.write(r.name + " " + r.frequency + " " +r.getLength() + " " + averrorsdist + " " + nrepl*r.frequency + " " + nins*r.frequency + " " + ndel*r.frequency + " ");
                    for (int i = 1; i <= homoplen; i++)
                        fw.write(errors[i]*r.frequency + " ");
                    fw.write(allerr*r.frequency + " " + nreplat*r.frequency + " " + nreplac*r.frequency + " " + nreplag*r.frequency + " " + nrepltc*r.frequency + " " + nrepltg*r.frequency + " " + nreplgc*r.frequency + "\n");
                    fw_un.write(r.name + " " + r.frequency + " " +r.getLength() + " " + averrorsdist + " " + nrepl + " " + nins + " " + ndel + " ");
                    for (int i = 1; i <= homoplen; i++)
                        fw_un.write(errors[i] + " ");
                    fw_un.write(allerr + " " + nreplat+ " " + nreplac+ " " + nreplag+ " " + nrepltc+ " " + nrepltg+ " " + nreplgc + "\n");
                        
            }
            DynamicOut.finishSteps();
            fw.close();
            fw_un.close();
            
        }
        void setFindErrorsSeglen(int i)
        {
            finderrorsseglen = i;
        }
        void fixDirection(double gapop, double gapext) throws IOException
        {
            String s_ref = "";
            for (Read r : reads)
                if (r.getLength() > s_ref.length())
                    s_ref = r.getNucl();
            for (Read r : reads)
            {
                System.out.println("Fixing direction: " + r.name);
                int bestmatch = 0;
                int bestdirect = -1;
                String nucl = r.getNucl();
                for (int i = 0; i <= 1; i++)
                {
                    Runtime run=Runtime.getRuntime();
                    Process p=null;
                    FileWriter fw_alin = new FileWriter("allign_input.fas");
                    fw_alin.write(">" + r.name + "\n");
                    fw_alin.write(nucl + "\n");
                    fw_alin.write(">reference" + "\n");
                    fw_alin.write(s_ref + "\n");
                    fw_alin.close();
                    String param = "ClustalW2//clustalw2 -INFILE=" +"allign_input.fas" + " -OUTFILE=" + "allign_output.fas";
                    param+= " -OUTPUT=FASTA -DNAMATRIX=IUB -GAPOPEN=" + gapop +" -GAPEXT=" + gapext + " -TYPE=DNA -PWDNAMATRIX=IUB -PWGAPOPEN=" + gapop+ " -PWGAPEXT=" + gapext;

                    p = run.exec(param);
                    BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
                    String s1;
                    while ((s1 = stdInput.readLine()) != null) {
    //                  System.out.println(s1);
                    }
                    if (p.exitValue()!= 0)
                    {
                        System.out.println("Error in allgnment program");
                    }
                    DataSet allignment = new DataSet("allign_output.fas",'c');
                    String al_read = "";
                    String al_ref = "";
                    for (Read r1 : allignment.reads)
                        if (!r1.name.contains("reference"))
                            al_read = r1.getNucl();
                        else
                            al_ref = r1.getNucl();
                    int match = 0;
                    for (int j = 0; j < al_read.length(); j++)
                        if (al_read.charAt(j) == al_ref.charAt(j))
                            match++;
                    if (match > bestmatch)
                    {
                        bestmatch = match;
                        bestdirect = i;
                    }
                    File f = new File("allign_input.fas");
                    f.delete();
                    f = new File("allign_output.fas");
                    f.delete();
                    f = new File("allign_input.dnd");
                    f.delete();
                    nucl = r.RevComp().getNucl();
                }
                if (bestdirect == 1)
                    r.nucl = r.RevComp().getNucl();
            }
        }

        public void PrintReads(String outfile) throws IOException
	{
		FileWriter fw = new FileWriter(outfile);
		for(Read r : reads)
                {
				for (int i = 0; i < r.frequency; i++)
					fw.write(">" + r.name + "_" + i + "\n" + r.getNucl() + "\n");
		}
		fw.close();
	}
        void fixDirection(String refFile, double gapop, double gapext) throws IOException
        {
            DataSet refds = new DataSet(refFile);
            String s_ref = refds.reads.get(0).getNucl();
            for (Read r : reads)
            {
                System.out.println("Fixing direction: " + r.name + "/" + reads.size());
                int bestmatch = 0;
                int bestdirect = -1;
                String nucl = r.getNucl();
                for (int i = 0; i <= 1; i++)
                {
                    Runtime run=Runtime.getRuntime();
                    Process p=null;
                    FileWriter fw_alin = new FileWriter("allign_input.fas");
                    fw_alin.write(">" + r.name + "\n");
                    fw_alin.write(nucl + "\n");
                    fw_alin.write(">reference" + "\n");
                    fw_alin.write(s_ref + "\n");
                    fw_alin.close();
                    String param = "ClustalW2//clustalw2 -INFILE=" +"allign_input.fas" + " -OUTFILE=" + "allign_output.fas";
                    param+= " -OUTPUT=FASTA -DNAMATRIX=IUB -GAPOPEN=" + gapop +" -GAPEXT=" + gapext + " -TYPE=DNA -PWDNAMATRIX=IUB -PWGAPOPEN=" + gapop+ " -PWGAPEXT=" + gapext;

                    p = run.exec(param);
                    BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
                    String s1;
                    while ((s1 = stdInput.readLine()) != null) {
    //                  System.out.println(s1);
                    }
                    if (p.exitValue()!= 0)
                    {
                        System.out.println("Error in allgnment program");
                    }
                    DataSet allignment = new DataSet("allign_output.fas",'c');
                    String al_read = "";
                    String al_ref = "";
                    for (Read r1 : allignment.reads)
                        if (!r1.name.contains("reference"))
                            al_read = r1.getNucl();
                        else
                            al_ref = r1.getNucl();
                    int match = 0;
                    for (int j = 0; j < al_read.length(); j++)
                        if (al_read.charAt(j) == al_ref.charAt(j))
                            match++;
                    if (match > bestmatch)
                    {
                        bestmatch = match;
                        bestdirect = i;
                    }
                    File f = new File("allign_input.fas");
                    f.delete();
                    f = new File("allign_output.fas");
                    f.delete();
                    f = new File("allign_input.dnd");
                    f.delete();
                    nucl = r.RevComp().getNucl();
                }
                if (bestdirect == 1)
                    r.nucl = r.RevComp().getNucl();
            }
        }
        void PrintReadsSeparateFiles() throws IOException
        {
            int i = 1;
            for (Read r : reads)
            {
                FileWriter fw = new FileWriter(this.file_name + "_Read" + i + ".fas");
                fw.write(">" + r.name + "\n" + r.nucl + "\n");
                fw.close();
                i++;
            }
            
        }
        void PrintReadsSeparateFilesKeepName(String addr) throws FileNotFoundException, IOException
        // <editor-fold defaultstate="collapsed" desc=" DESCRIPTION ">
    {
        FileReader fr = new FileReader(addr);
        BufferedReader br = new BufferedReader(fr);
        int i = 1;
        String s = br.readLine();
        String nucl = new String();
        String name = s;
        s = br.readLine();
        while (s != null) {
            if (s.length() == 0) {
                s = br.readLine();
                continue;
            }
            if (s.charAt(0) == '>') {
                FileWriter fw = new FileWriter(addr + "_Read" + i + ".fas");
                fw.write(name + "\n" + nucl + "\n");      
                fw.close();
                nucl = "";
                name = s;
                i++;
            } else {
                nucl += s.toUpperCase();
            }
            s = br.readLine();
        }
        FileWriter fw = new FileWriter(addr + "_Read" + i + ".fas");
        fw.write(name + "\n" + nucl + "\n");
        fw.close();
        br.close();
        fr.close();
        
    }
// </editor-fold>
        void printMostFreqReads(String addr, int n) throws IOException
        {
            Collections.sort(reads, new ReadFreqComparator());
            FileWriter fw = new FileWriter(addr);
            for (int i = 1; i <= n; i++)
            {
                fw.write(">Read" + i + "_freq" + reads.get(i).frequency + "\n" + reads.get(i).getNucl() + "\n");
            }
            fw.close();
        }
        Read findMostFreqRead()
                // <editor-fold defaultstate="collapsed" desc=" DESCRIPTION ">
        {
            Read maxr = new Read();
            double maxfreq = 0;
            for (Read r : this.reads)
                if (r.frequency > maxfreq)
                {
                    maxfreq = r.frequency;
                    maxr = r;
                }
            return maxr;
        }
        // </editor-fold>
        void printpairdistReads(String addr) throws IOException
                // <editor-fold defaultstate="collapsed" desc=" DESCRIPTION ">
        {
            FileWriter pairdist = new FileWriter(addr);
            int[][] pd = new int[reads.size()][reads.size()];

            double gapop = 6;
            double gapext = 6;
            for (int i = 0; i < reads.size(); i++)
                for (int j = i+1; j < reads.size(); j++)
                {
                    DynamicOut.printStep("Calculating pair (" + i + "," + j + ") //" + reads.size());
                    Runtime run=Runtime.getRuntime();
                    Process p=null;
                    FileWriter fw_alin = new FileWriter("allign_input.fas");
                    fw_alin.write(">first" + "\n");
                    fw_alin.write(reads.get(i).getNucl() + "\n");
                    fw_alin.write(">second" + "\n");
                    fw_alin.write(reads.get(j).getNucl() + "\n");
                    fw_alin.close();
                    String param = "ClustalW2//clustalw2 -INFILE=" +"allign_input.fas" + " -OUTFILE=" + "allign_output.fas";
                    param+= " -OUTPUT=FASTA -DNAMATRIX=IUB -GAPOPEN=" + gapop +" -GAPEXT=" + gapext + " -TYPE=DNA -PWDNAMATRIX=IUB -PWGAPOPEN=" + gapop+ " -PWGAPEXT=" + gapext;

                    p = run.exec(param);
                    BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
                    String s1;
                    while ((s1 = stdInput.readLine()) != null) {
        //              System.out.println(s1);
                    }
                    if (p.exitValue()!= 0)
                    {
                        System.out.println("Error in allgnment program");
                    }
                    DataSet allignment = new DataSet("allign_output.fas",'c');
                    int l = allignment.reads.get(0).getLength();
                    int mismatch = 0;
                    for (int k = 0; k < l; k++)
                        if (allignment.reads.get(0).getNucl().charAt(k) != allignment.reads.get(1).getNucl().charAt(k))
                            mismatch++;
                    pd[i][j] = mismatch;
                    pd[j][i] = mismatch;
                    File f = new File("allign_input.fas");
                    f.delete();
                    f = new File("allign_output.fas");
                    f.delete();
                    f = new File("allign_input.dnd");
                    f.delete();
                }
            DynamicOut.finishSteps();
            for (int i = 0; i < reads.size(); i++)
            {
                for (int j = 0; j < reads.size(); j++)
                    pairdist.write(pd[i][j] + " ");
                pairdist.write("\n");
            }
            pairdist.close();
        }
        // </editor-fold>
}
