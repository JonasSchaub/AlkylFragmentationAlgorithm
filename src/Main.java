
public class Main {

	public static void main(String[] args) {
		String molecule = "CC(CCC)(CC)CCCCC(CC)C(CCC)CCCCCCCC(C)(C)C";
		Fragmenter frag = new Fragmenter(molecule);
		frag.printInformation(3);
	}

}
