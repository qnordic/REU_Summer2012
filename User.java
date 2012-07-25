import java.util.Set;
import java.util.TreeSet;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;

public class User {
	public Set<Scite> paper ;
	public Set<Date> activeDates ;
	
	public User(String paper, String date) {
		
	}
	
	public User(Scite paper, String date) {
		this.paper = paper ;
		activeDates = new TreeSet<Date>() ;
		activeDates.add(parseDate(date)) ;
	}
	
	public void addDate(String date) {
		activeDates.add(parseDate(date)) ;
	}
	
	private Date parseDate(String input) {
		Date date = null ;
		DateFormat formatter ;
		try {
			formatter = new SimpleDateFormat("MM/dd/yyyy HH:mm") ;
			date = (Date) formatter.parse(input) ;
		} catch (ParseException e) {
			System.out.println("Exception: " + e) ;
		}
		return date ;
	}
}
