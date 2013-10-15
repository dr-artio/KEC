/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package errorcorrection;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 *
 * @author kki8
 */
public class ForDavid {
     public static void main(String[] args) throws FileNotFoundException, IOException
     {
                FileReader fr = new FileReader("2.GAC.454Reads.fna");
		BufferedReader br = new BufferedReader(fr);

                int num = 0;
                int len = 0;

		String s = br.readLine();
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
				//	reads.add(new Read(nucl, name));

//					System.out.println(name);
                                        len+=nucl.length();
                                        num++;
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

                len+=nucl.length();
                num++;
		fr.close();
                System.out.println((double) len/num);
                System.out.println(num);
     }

}
