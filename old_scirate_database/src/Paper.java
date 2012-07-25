import java.util.* ;

public class Paper implements Comparable<Paper> {
	public String paperID ;
	public Date firstScite ;
	
	public Paper(String paperID, Date firstScite) {
		
		this.paperID = paperID ;
		this.firstScite = firstScite ;
	}
	
	public int compareTo(Paper other) {
		if (other.firstScite.before(this.firstScite))
			return 1 ;
		else if (other.firstScite.equals(this.firstScite))
			return 0 ;
		else
			return -1 ;
	}
}