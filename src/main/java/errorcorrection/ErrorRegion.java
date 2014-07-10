package ErrorCorrection;

public class ErrorRegion extends Sequence {
	int begin;
	int end;
	int length;
	ErrorRegion(String s, int i, int j, int l)
	{
		super(s);
		begin = i;
		end = j;
		length = l;
	}
}
