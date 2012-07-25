import java.io.* ;
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
		int index = -1 ;
		
		Map<String,ArrayList<Integer>> userPaper = new TreeMap<String,ArrayList<Integer>>() ;
		Map<String,Integer> papers = new HashMap<String,Integer>() ;
		
		// Creates the matrix for use without sorting by date.
		while (fileScan.hasNextLine()) {
			tempLine = fileScan.nextLine().split(",") ;
			if (!tempLine[3].contains("1900-01-01")) {
				if (!papers.containsKey(tempLine[1]))
					papers.put(tempLine[1], index++) ;
				if (!userPaper.containsKey(tempLine[2]))
					userPaper.put(tempLine[2], new ArrayList<Integer>()) ;
				userPaper.get(tempLine[2]).add(papers.get(tempLine[1])) ;
			}
		}
		
		Set<String> users = userPaper.keySet() ;
		System.out.println("Finished sorting values") ;
		
		try {
			FileWriter fstream = new FileWriter("oldDataMatrix.txt");
			BufferedWriter out = new BufferedWriter(fstream);
			for (Map.Entry<String, ArrayList<Integer>> entry : userPaper.entrySet())
				for (int i = 0; i <= papers.size(); i++)
					if (i == papers.size())
						out.write("\n") ;
					else if (entry.getValue().contains(i))
						out.write("1,") ;
					else
						out.write("0,") ;
			out.close();
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
		}
		
		try { 
			FileWriter fstream = new FileWriter("oldDataUsers.txt");
			BufferedWriter out = new BufferedWriter(fstream);
			for (String temp : users)
				out.write(temp + "\n") ;
			out.close();
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
		}
		try { 
			FileWriter fstream = new FileWriter("oldDataPapers.txt");
			BufferedWriter out = new BufferedWriter(fstream);
			for (Map.Entry<String,Integer> entry : papers.entrySet())
				out.write(entry.getKey() + ", " + entry.getValue() + "\n") ;
			out.close();
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
		}
		System.out.println("Finished writing files") ;
	}
}