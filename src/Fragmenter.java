import java.util.ArrayList;
import java.util.List;

public class Fragmenter {
	private String mol;
	private List<Integer> secondaryCarbon;
	private List<Integer> tertiaryCarbon;
	private List<Fragment> extractedFragments;
	
	public Fragmenter(String molecule) {
		mol = molecule;
		secondaryCarbon = branches(false);
		tertiaryCarbon = branches(true);
		removeBranches();
	}
	
	private List<Integer> branches(boolean tert) {
		List<Integer> branchIndex = new ArrayList<>();
		List<Integer> bracketIndex = new ArrayList<>();
		for (int i=0; i<mol.length(); i++) {
			if (mol.charAt(i) == '(') {
				bracketIndex.add(i);
			}
			if (mol.charAt(i) == ')') {
				if (tert && mol.charAt(i+1) == '(') {
					branchIndex.add(bracketIndex.get(bracketIndex.size()-1)-1);
				}
				if (!tert && mol.charAt(i+1) != '(' && mol.charAt(bracketIndex.get(bracketIndex.size()-1)-1) != ')') {
					branchIndex.add(bracketIndex.get(bracketIndex.size()-1)-1);
				}
				bracketIndex.remove(bracketIndex.size()-1);
			}
		}
		return branchIndex;
	}
	
	private int countCarbons() {
		int n = 0;
		for (int i=0; i<mol.length(); i++) {
			if (mol.charAt(i) == 'C') {
				n++;
			}
		}
		return n;
	}
	
	private void removeBranches() {
		extractedFragments = new ArrayList<>();
		int start = 0;
		for (int i=0; i<mol.length(); i++) {
			if (secondaryCarbon.contains(i) || tertiaryCarbon.contains(i) || mol.charAt(i) == ')' || i == mol.length()-1) {
				
				//+(int)(i/(mol.length()-1)) → becomes 1 at the end of String, so that last char will be included in substring 
				String newFragment = mol.substring(start, i+(int)(i/(mol.length()-1)));
				if (newFragment != "" && containsFragment(extractedFragments, newFragment)) {
					extractedFragments.get(indexOfFragment(extractedFragments, newFragment)).addIndex(start);
				}
				if (newFragment != "" && !containsFragment(extractedFragments, newFragment)) {
					extractedFragments.add(new Fragment(newFragment, start));
				}
			}
			if (mol.charAt(i) == '(' || mol.charAt(i) == ')' && mol.charAt(i+1) != '(') {
				start = i+1;
			}
		}
	}

	private List<Fragment> splitFragments(int alkSize) {
		List<Fragment> splitFrags = new ArrayList<>();
		for (int i=0; i<extractedFragments.size(); i++) {
			if (extractedFragments.get(i).getFragment().length() >= alkSize) {
				if (!containsFragment(splitFrags, "C".repeat(alkSize))) {
					splitFrags.add(new Fragment("C".repeat(alkSize)));
				}
				if (containsFragment(splitFrags, "C".repeat(alkSize))) {
					for (int j=0; j<(int)(extractedFragments.get(i).getFragment().length()/alkSize); j++) {
						for (int k=0; k<extractedFragments.get(i).getIndices().size(); k++) {
							splitFrags.get(indexOfFragment(splitFrags, "C".repeat(alkSize))).addIndex(extractedFragments.get(i).getIndices().get(k)+j*alkSize);
						}
					}
				}
			}
			if (!containsFragment(splitFrags, "C".repeat(extractedFragments.get(i).getFragment().length()%alkSize)) && extractedFragments.get(i).getFragment().length()%alkSize != 0) {
				splitFrags.add(new Fragment("C".repeat(extractedFragments.get(i).getFragment().length())));
			}
			if (containsFragment(splitFrags, "C".repeat(extractedFragments.get(i).getFragment().length()%alkSize)) && extractedFragments.get(i).getFragment().length()%alkSize != 0) {
				for (int k=0; k<extractedFragments.get(i).getIndices().size(); k++) {
					splitFrags.get(indexOfFragment(splitFrags, "C".repeat(extractedFragments.get(i).getFragment().length()%alkSize))).addIndex(extractedFragments.get(i).getIndices().get(k)+alkSize*(int)(extractedFragments.get(i).getFragment().length()/alkSize));
				}	
			}
		}
		return splitFrags;
	}
	
	private boolean containsFragment(List<Fragment> l, String frag) {
		boolean contained = false;
		for (Fragment li : l) {
			if (li.getFragment().equals(frag)) {
				contained = true;
			}
		}
		return contained;
	}
	
	private int indexOfFragment(List<Fragment> l, String frag) {
		int index = -1;
		for (int i=0; i<l.size(); i++) {
			if (l.get(i).getFragment().equals(frag)) {
				index = i;
			}
		}
		return index;
	}
	
	private void listPrinter(List<Integer> l) {
		for (int i=0; i<l.size(); i++) {
			System.out.print(l.get(i));
			if (i < l.size()-1) {
				System.out.print(", ");
			}
		}
	}
	
	private void listPrinterFrag(List<Fragment> l) {
		for (int i=0; i<l.size(); i++) {
			System.out.print(l.get(i).getFragment() + "(");
			for (int j=0; j<l.get(i).getIndices().size(); j++) {
				System.out.print(l.get(i).getIndices().get(j));
				if (j < l.get(i).getIndices().size()-1) {
					System.out.print(", ");
				}
			}
			if (i < l.size()-1) {
				System.out.print("), ");
			}
		}
		System.out.print(')');
	}
	
	public void printInformation(int alkSize) {
		System.out.print("Molecule: " + mol);
		System.out.print("\n\nnumber of carbon atoms: " + Integer.toString(countCarbons()));
		System.out.print("\n\nIndices of ...\n\t... secondary carbons: ");
		listPrinter(secondaryCarbon);		
		System.out.print("\n\t... tertiary carbons: ");
		listPrinter(tertiaryCarbon);
		System.out.print("\n\nextracted fragments: ");
		listPrinterFrag(extractedFragments);
		System.out.print(String.format("\n\nsplit fragments (%sC): ",alkSize));
		listPrinterFrag(splitFragments(alkSize));
	}
}
