import java.util.* ;
import java.text.* ;
public class Paper implements Comparable<Paper> {
	public String paperID ;
	public Date firstScite ;
	
	public Paper(String paperID, String firstScite) {
		this.paperID = paperID ;
		DateFormat formatter ;
		try {
			formatter = new SimpleDateFormat("MM/dd/yyyy HH:mm") ;
			this.firstScite = (Date) formatter.parse(firstScite) ;
		} catch (ParseException e) {
			System.out.println("Exception: " + e) ;
		}
	}
	
	public int compareTo(Paper other) {
		return (int) (this.firstScite.getTime() - other.firstScite.getTime()) ;
	}
}