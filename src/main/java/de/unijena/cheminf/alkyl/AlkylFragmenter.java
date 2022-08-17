package de.unijena.cheminf.alkyl;
/*
 * MIT License
 *
 * Copyright (c) 2022 Andreas Freitag, Jonas Schaub, Felix Bänsch, Achim Zielesny, Christoph Steinbeck
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

import org.openscience.cdk.interfaces.*;

import java.util.ArrayList;
import java.util.List;

public class AlkylFragmenter {
    private IAtomContainer mol;

    public AlkylFragmenter(){}

    public void setMolecule(IAtomContainer m) {
        mol = m;
    }
    public void startAlkylFragmentation() {
        List<List<Integer>> bc = branchCutter(mol);
        System.out.println(bc);
        System.out.println(chainCutter(mol, bc, 5, 0, false));
    }

    private List<IAtomContainer> genAtomContainer(List<List<Integer>> indicesList) throws CloneNotSupportedException {
        List<IAtomContainer> molList = new ArrayList<>();
        for (List<Integer> il : indicesList) {
            IAtomContainer molecule = mol.clone();
            for (int i=molecule.getAtomCount()-1; i<=0; i++) {
                if (!il.contains(i)) {
                    molecule.removeAtom(i);
                }
            }
            molList.add(molecule.clone());
        }
        return molList;
    }

    /**
     * This method separates the shorter branches from the longer ones and returns a list of linear molecule fragments
     * in shape of lists of atom indices.
     * @param molecule
     * @return
     */
    public List<List<Integer>> branchCutter (IAtomContainer molecule) {
        List<List<Integer>> chains = new ArrayList<>();
        List<List<Integer>> chainList = new ArrayList<>();
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
                chainList.add(new ArrayList<>());
                chainList.get(chainList.size()-1).add(i);
            }
        }
        boolean searching = true;
        while (searching) {
            int chainListIndex = 0;
            while (chainListIndex < chainList.size()) {
                List<Integer> chainListAtIndex = chainList.get(chainListIndex);
                if (connections[chainListAtIndex.get(chainListAtIndex.size() - 1)].size() > 0 &&
                        chainList.size() > 1) {
                    chainListAtIndex.add(connections[chainListAtIndex.get(chainListAtIndex.size() - 1)].get(0));
                    connections[chainListAtIndex.get(chainListAtIndex.size() - 1)]
                            .remove(chainListAtIndex.get(chainListAtIndex.size() - 2));
                    if (connections[chainListAtIndex.get(chainListAtIndex.size() - 1)].size() > 1) {
                        chains.add(reverseList(chainListAtIndex));
                        chainList.remove(chainListAtIndex);
                        chainListIndex--;
                    }
                } else {
                    if (chainList.size() == 2) {
                        chains.add(chainList.get(0));
                        chainList.remove(0);
                    }
                    chains.get(chains.size() - 1).remove(chains.get(chains.size() -1).size() - 1);
                    chains.get(chains.size() - 1).addAll(reverseList(chainList.get(0)));
                    chainList.remove(0);
                    chainListIndex--;
                    searching = false;
                }
                chainListIndex++;
            }
        }
        return chains;
    }

    private List<Integer> reverseList(List<Integer> list) {
        List<Integer> revList = new ArrayList<>();
        for (int li : list) {
            revList.add(0,li);
        }
        return revList;
    }

    /**
     * This method cuts linear molecule fragments into chains of equal length.
     * @param molecule
     * @param chains
     * @param minCut
     * @param maxCut
     * @param preserveCarbon
     * @return
     */
    public List<List<Integer>> chainCutter(IAtomContainer molecule, List<List<Integer>> chains,
                                           int minCut, int maxCut, boolean preserveCarbon) {
        List<List<Integer>> rest = new ArrayList<>();
        int chainsIndex = 0;
        while (chainsIndex < chains.size() && chains.size() > 0) {
            List<Integer> chainsAtIndex = chains.get(chainsIndex);
            int index = 1;
            int branchRest = 0;
            while (index-branchRest < 2 && index != chains.size()-1) {
                if (index == 1 && preserveCarbon) branchRest++;
                else if (molecule.getBond(molecule.getAtom(chainsAtIndex.get(branchRest)),
                        molecule.getAtom(chainsAtIndex.get(index))).getOrder() != IBond.Order.SINGLE) branchRest++;
                index++;
            }
            if (index-branchRest > 0) {
                if (chainsAtIndex.size() - branchRest > minCut) {
                    if (branchRest > 0) rest.add(chainsAtIndex.subList(0, branchRest));
                    chainsAtIndex = chainsAtIndex.subList(branchRest, chainsAtIndex.size());
                    if (chainsIndex < chains.size()-1) chainsAtIndex.remove(0);
                    int shift = 0;
                    index = 0;
                    if (maxCut > 0) {
                        while ((index+1) * maxCut + shift <= chainsAtIndex.size()) {
                            int shift0 = shift;
                            boolean reverseShift = false;
                            if ((index + 1) * maxCut + shift + 1 < chainsAtIndex.size()) {
                                while (molecule.getBond(molecule.getAtom(chainsAtIndex.get((index + 1) * maxCut + shift)),
                                        molecule.getAtom(chainsAtIndex.get((index + 1) * maxCut + shift + 1))).getOrder() != IBond.Order.SINGLE &&
                                        (index + 1) * maxCut + shift < chainsAtIndex.size()) {
                                    if (shift0 - shift < maxCut - minCut && !reverseShift) shift--;
                                    if (shift0 - shift >= maxCut - minCut && !reverseShift) {
                                        reverseShift = true;
                                        shift = shift0;
                                    }
                                    if (reverseShift) shift++;
                                }
                            }
                            chains.add(0, chainsAtIndex.subList(index*maxCut+shift, (index+1)*maxCut+shift));
                            chainsIndex++;
                            index++;
                        }
                        if (chainsAtIndex.size()/maxCut < (float) chainsAtIndex.size()/maxCut)
                            chains.add(0, chainsAtIndex.subList(index*maxCut+shift, chainsAtIndex.size()));
                        chains.remove(chainsAtIndex);

                    } else {
                        while ((index + 1) * minCut + shift <= chainsAtIndex.size()) {
                            if ((index + 1) * minCut + shift + 1 < chainsAtIndex.size()) {
                                while (molecule.getBond(molecule.getAtom(chainsAtIndex.get((index + 1) * minCut + shift)),
                                        molecule.getAtom(chainsAtIndex.get((index + 1) * minCut + shift + 1))).getOrder() != IBond.Order.SINGLE &&
                                        (index + 1) * minCut + shift < chainsAtIndex.size()) {
                                    shift++;
                                }
                            }
                            chains.add(0, chainsAtIndex.subList(index * minCut + shift, (index + 1) * minCut + shift));
                            chainsIndex++;
                            index++;
                        }
                        if (chainsAtIndex.size()/minCut < (float) chainsAtIndex.size()/minCut)
                            rest.add(chainsAtIndex.subList(index*minCut+shift-1, chainsAtIndex.size()));
                        chains.remove(chainsAtIndex);
                        chainsIndex--;
                    }
                } else {
                    rest.add(chainsAtIndex);
                    chains.remove(chainsAtIndex);
                    chainsIndex--;
                }
            }
            chainsIndex++;
        }
        while (rest.size() > 0) {
            chainsIndex = 0;
            while (chainsIndex < chains.size()) {
                List<Integer> chainsAtIndex = chains.get(chainsIndex);
                int restIndex = 0;
                boolean isCombined = false;
                while (restIndex < rest.size() && !isCombined) {
                    List<Integer> restAtIndex = rest.get(restIndex);
                    if (chainsAtIndex.contains(restAtIndex.get(0))) {
                        List<Integer> combine = new ArrayList<>(chainsAtIndex);
                        combine.addAll(restAtIndex.subList(1, restAtIndex.size()));
                        chains.add(combine);
                        chains.remove(chainsAtIndex);
                        rest.remove(restAtIndex);
                        restIndex++;
                        isCombined = true;
                    }
                    restIndex++;
                }
                if (!isCombined) chainsIndex++;
            }
        }
        return chains;
    }
}