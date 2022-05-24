package de.unijena.cheminf.alkyl;
/*
 * MIT License
 *
 * Copyright (c) 2022 Andreas Freitag, Jonas Schaub, Achim Zielesny, Christoph Steinbeck
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

public class AlkylFragmenter {
    private IAtomContainer mol;
    private List<List<Integer>> ccFragments;    //ccFragments - main and side branch fragments
    private List<List<Integer>> tqFragments;    //tqFragments - tertiary and quaternary and chain fragments


    public AlkylFragmenter (IAtomContainer m) throws CDKException {
        mol = m;
        System.out.println(cutChains(mol));

    }

    public List<Integer> longestChainIndices (IAtomContainer molecule) {
        return cutChains(molecule).get(0);
    }

    public List<List<Integer>> cutChains (IAtomContainer molecule) {
        /*
        This method determines the longest chain in a hydrocarbon molecule by following chain paths starting with the
        terminal carbons atom by atom all at the same time. When any path reaches a fork it will become a side
        chain and its list will be deleted. This way the two longest paths will never reach a fork, but will meet each
        other in the middle.
         */
        List<List<Integer>> chains = new ArrayList<>();

        List<List<Integer>> chainLists = new ArrayList<>();
        List<Integer>[] connections = new ArrayList[molecule.getAtomCount()];
        for (int i=0; i<connections.length; i++) {
            connections[i] = new ArrayList<>();
        }
        for (IBond bond : molecule.bonds()) {
            connections[bond.getAtom(0).getIndex()].add(bond.getAtom(1).getIndex());
            connections[bond.getAtom(1).getIndex()].add(bond.getAtom(0).getIndex());
        }
        for (int i=0; i<connections.length; i++) {
            if (connections[i].size() == 1) {
                chainLists.add(new ArrayList<>());
                chainLists.get(chainLists.size()-1).add(i);
            }
        }

        boolean searching = true;
        while (searching) {
            int cli = 0;
            while (cli < chainLists.size()) {
                List<Integer> cl = chainLists.get(cli);
                if (connections[cl.get(cl.size()-1)].size() > 0) {
                    cl.add(connections[cl.get(cl.size() - 1)].get(0));
                    connections[cl.get(cl.size() - 1)].remove((Integer) cl.get(cl.size() - 2));
                    if (connections[cl.get(cl.size() - 1)].size() == 0) {
                        if (Objects.equals(chains.get(0).get(chains.get(0).size() - 1), cl.get(cl.size() - 1))) {
                            searching = false;
                        }
                        cl.remove(cl.size() - 1);
                        chains.add(0, cl);
                        chainLists.remove(cl);
                        cli--;
                    }
                    if (connections[cl.get(cl.size() - 1)].size() > 1){
                        cl.remove(cl.size()-1);
                        chains.add(cl);
                        chainLists.remove(cl);
                        cli--;
                    }
                } else {
                    chains.get(0).addAll(cl);
                    chainLists.remove(cl);
                    cli--;
                    searching = false;
                }
                cli++;
            }
        }
        return chains;
    }

}