import java.io.* ;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.* ;

public class oldScirateData {
	public static void main(String[] args) {
		
		File data = new File("votes.csv") ;
		Scanner fileScan = null ;
		
		try {
			fileScan = new Scanner(data) ;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		String[] tempLine = null ;
		
		Map<String,Date> firstPaperScite = new HashMap<String,Date>() ;
		Map<String,User> users = new HashMap<String,User>() ;
		User tempUser = null ;
		
		while (fileScan.hasNextLine()) {
			tempLine = fileScan.nextLine().split(",") ;
			if (!tempLine[3].contains("1900-01-01")) {
				if (!firstPaperScite.containsKey(tempLine[1]))
					firstPaperScite.put(tempLine[1],parseDate(tempLine[3])) ;
				else
					if (firstPaperScite.get(tempLine[1]).after(parseDate(tempLine[3])))
						firstPaperScite.put(tempLine[1],parseDate(tempLine[3])) ;
				
				tempUser = new User(tempLine[2], tempLine[1], tempLine[3]) ;
				if (!users.containsKey(tempLine[2]))
					users.put(tempLine[2], tempUser) ;
				else
					users.get(tempLine[2]).addScite(tempLine[1], tempLine[3]) ;
			}
		}
		
		PriorityQueue<Paper> papers = new PriorityQueue<Paper>() ;
		for (Map.Entry<String,Date> entry : firstPaperScite.entrySet())
			papers.add(new Paper(entry.getKey(), entry.getValue())) ;
		
		ArrayList<Paper> paperIndex = new ArrayList<Paper>() ;
		while (!papers.isEmpty())
			paperIndex.add(papers.remove()) ;
		
		System.out.println("Finished sorting values") ;
		
		printFiles("matrix", users, paperIndex) ;
		printFiles("users", users, null) ;
		printFiles("papers", null, paperIndex) ;
		
		System.out.println("Finished writing files") ;
	}
	
	public static Date parseDate(String input) {
		Date date = null ;
		DateFormat formatter ;
		try {
			formatter = new SimpleDateFormat("MM/dd/yyyy") ;
			date = (Date) formatter.parse(input) ;
		} catch (ParseException e) {
			System.out.println("Exception: " + e) ;
		}
		return date ;
	}
	
	public static void printFiles(String option, Map<String,User> userPaper, ArrayList<Paper> papers) {
		
		if (option.equalsIgnoreCase("matrix"))
			try {
				FileWriter fstream = new FileWriter("oldDataMatrix.txt") ;
				BufferedWriter out = new BufferedWriter(fstream) ;
				
				for (Map.Entry<String,User> entry : userPaper.entrySet()) {
					for (int i = 0; i <= papers.size(); i++)
						if (i == papers.size())
							out.write("\n") ;
						else if (entry.getValue().paperIDs.contains(papers.get(i).paperID))
							out.write("1,") ;
						else if (entry.getValue().activeDays.contains(papers.get(i).firstScite))
							out.write("-1,") ;
						else
							out.write("0,") ;
				}
				out.close() ;
				
			} catch (Exception e) {
				System.err.println("Error: " + e.getMessage()) ;
			}
		
		if (option.equalsIgnoreCase("users"))
			try { 
				FileWriter fstream = new FileWriter("oldDataUsers.txt") ;
				BufferedWriter out = new BufferedWriter(fstream) ;
				
				for (Map.Entry<String,User> entry : userPaper.entrySet())
					out.write(entry.getValue().userID + "\n") ;
				out.close() ;
				
			} catch (Exception e) {
				System.err.println("Error: " + e.getMessage()) ;
			}
		
		if (option.equalsIgnoreCase("papers"))
			try { 
				FileWriter fstream = new FileWriter("oldDataPapers.txt") ;
				BufferedWriter out = new BufferedWriter(fstream) ;
				
				for (int i = 0; i < papers.size(); i++)
					out.write(papers.get(i).paperID + ", " + papers.get(i).firstScite + "\n") ;
				out.close() ;
				
			} catch (Exception e) {
				System.err.println("Error: " + e.getMessage()) ;
			}
	}
}