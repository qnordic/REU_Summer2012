import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;


public class Scite {
	public String paper ;
	public Date date ;
	
	public Scite (String paper, String date) {
		this.paper = paper ;
		DateFormat formatter ;
		try {
			formatter = new SimpleDateFormat("MM/dd/yyyy HH:mm") ;
			this.date = (Date) formatter.parse(date) ;
		} catch (ParseException e) {
			System.out.println("Exception: " + e) ;
		}
	}
}
