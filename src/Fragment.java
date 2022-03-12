import java.util.ArrayList;
import java.util.List;

public class Fragment {
	private String frag;
	private List<Integer> indices;
	
	public Fragment(String f) {
		frag = f;
		indices = new ArrayList<>();
	}
	public Fragment(String f, int i) {
		frag = f;
		indices = new ArrayList<>();
		indices.add(i);
	}
	public void addIndex(int i) {
		indices.add(i);
	}
	public String getFragment() {
		return frag;
	}
	public List<Integer> getIndices() {
		return indices;
	}
}
