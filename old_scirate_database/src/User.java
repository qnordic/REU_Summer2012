import java.util.Set;
import java.util.TreeSet;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;

public class User {
	public String userID ;
	public Set<Date> activeDays ;
	public Set<String> paperIDs ;
	
	public User(String ID, String paperID, String date) {
		this.userID = ID ;
		paperIDs = new TreeSet<String>() ;
		paperIDs.add(paperID) ;
		activeDays = new TreeSet<Date>() ;
		activeDays.add(parseDate(date)) ;
	}
	
	public void addScite(String paperID, String date) {
		paperIDs.add(paperID) ;
		activeDays.add(parseDate(date)) ;
	}
	
	private Date parseDate(String input) {
		Date date = null ;
		DateFormat formatter ;
		if (input.contains("1900-01-01"))
			try {
				formatter = new SimpleDateFormat("yyyy-MM-dd") ;
				date = (Date) formatter.parse(input) ;
			} catch (ParseException e) {
				System.out.println("Exception: " + e) ;
			}
		else
			try {
				formatter = new SimpleDateFormat("MM/dd/yyyy") ;
				date = (Date) formatter.parse(input) ;
			} catch (ParseException e) {
				System.out.println("Exception: " + e) ;
			}
		return date ;
	}
}