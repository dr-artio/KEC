package ErrorCorrection;


public class Calculator {

	public int editDistance(String s, String t)
	{
		int m[][] = new int[s.length()+1][t.length()+1];
		for (int i = 0; i<=s.length();i++)
			m[i][0] = 0;
		for (int j = 0; j<=t.length();j++)
			m[0][j] = 0;
		for (int i = 1; i <= s.length(); i++)
			for (int j = 0; j <= t.length(); j++)
			{
				int a,b,c;
				if (s.charAt(i) == t.charAt(j))
					a = m[i-1][j-1] + 1;
				else
					a = 0;
				b = m[i-1][j];
				c = m[i][j-1];
				if (a > b)
					m[i][j] = Math.max(a, c);
				else
					m[i][i] = Math.max(b, c);				
			}
		return m[s.length()][t.length()];
	}
}
