package ErrorCorrection;

public class Correction {
	
//	D-deletion, R - replacement, I - insertion after
	char type;
	String replacement;
	Correction(char c, String s)
	{
		type = c;
		replacement = s;
	}
	Correction(char c, char s)
	{
		type = c;
		replacement =""+ s;
	}
	Correction(char c)
	{
		type = c;
	}

}
