/*
 * MIT License
 *
 * Copyright (c) 2022 Andreas Freitag, Jonas Schaub, Felix Baensch, Achim Zielesny, Christoph Steinbeck
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

package de.unijena.cheminf.alkyl;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.ArrayList;
import java.util.List;

public class AlkylFragmenter {
    private IAtomContainer molecule;
    private List<Integer>[] connections;
    private List<List<Integer>> branches;
    private List<List<Integer>> fragmentsIndices;
    private List<IAtomContainer> fragmentsAtomContainer;
    private List<List<Integer>> rest;
    private int minCut;
    private int maxCut;
    private boolean isPreservingTertiaryQuaternaryCarbons;

    public AlkylFragmenter(){}
    public void setMolecule(IAtomContainer aMolecule) {
        this.molecule = aMolecule;
    }
    public IAtomContainer getMolecule () {
        return this.molecule;
    }
    public void fragmentationSettings(int aMinCut, int aMaxCut, boolean aIsPreservingTertiaryQuaternaryCarbons) throws CloneNotSupportedException, CDKException {
        this.minCut = aMinCut;
        this.maxCut = aMaxCut;
        this.isPreservingTertiaryQuaternaryCarbons = aIsPreservingTertiaryQuaternaryCarbons;
        cutBranches();
        cutChains();
        System.out.println(this.branches);
        System.out.println(this.fragmentsIndices);
        System.out.println(this.rest);
        makeCorrections();
        genAtomContainer(this.fragmentsIndices);
        System.out.println(this.fragmentsIndices);
    }

    private void genAtomContainer(List<List<Integer>> anIndicesList) throws CloneNotSupportedException, CDKException {
        this.fragmentsAtomContainer = new ArrayList<>(anIndicesList.size());
        for (List<Integer> tmpListItem : anIndicesList) {
            IAtomContainer tmpMoleculeFragment = new AtomContainer();
            for (IBond tmpBond : molecule.bonds()) {
                if (tmpListItem.contains(tmpBond.getAtom(0).getIndex()) && tmpListItem.contains(tmpBond.getAtom(1).getIndex())) {
                    tmpMoleculeFragment.addBond(tmpBond);
                }
            }
            for (int tmpAtomIndex : tmpListItem) {
                tmpMoleculeFragment.addAtom(this.molecule.getAtom(tmpAtomIndex));
            }
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMoleculeFragment);
            CDKHydrogenAdder.getInstance(tmpMoleculeFragment.getBuilder()).addImplicitHydrogens(tmpMoleculeFragment);
            this.fragmentsAtomContainer.add(tmpMoleculeFragment);
        }
    }

    private IAtomContainer saturateMolecule (IAtomContainer aMolecule) {
        for (int i=0; i<aMolecule.getAtomCount(); i++) {
            for (int j=0; j<aMolecule.getAtom(i).getImplicitHydrogenCount(); j++) {
                IAtom hydrogen = aMolecule.getAtom(i).getBuilder().newInstance(IAtom.class, "H");
                hydrogen.setAtomTypeName("H");
                hydrogen.setImplicitHydrogenCount(0);
                aMolecule.addAtom(hydrogen);
                aMolecule.addBond(aMolecule.getAtom(i).getBuilder().newInstance(IBond.class, aMolecule.getAtom(i), hydrogen, IBond.Order.SINGLE));
            }
            aMolecule.getAtom(i).setImplicitHydrogenCount(0);
        }
        return aMolecule;
    }

    public List<IAtomContainer> getIAtomContainer () {
        return this.fragmentsAtomContainer;
    }

    public List<String> getSmiles () throws CDKException {
        SmilesGenerator tmpSmilesGenerator = SmilesGenerator.generic();
        List<String> tmpSmilesList = new ArrayList<>(this.fragmentsAtomContainer.size());
        for (IAtomContainer tmpFragment : this.fragmentsAtomContainer) {
            tmpSmilesList.add(tmpSmilesGenerator.create(tmpFragment));
        }
        return tmpSmilesList;
    }

    /**
     * This method dissects the molecule into individual branches.
     */
    public void cutBranches () {
        this.branches = new ArrayList<>();
        List<List<Integer>> tmpChainList = new ArrayList<>();
        this.connections = new ArrayList[this.molecule.getAtomCount()];
        for (int i=0; i<this.connections.length; i++) {
            this.connections[i] = new ArrayList<>();
        }
        for (IBond bond : this.molecule.bonds()) {
            this.connections[bond.getAtom(0).getIndex()].add(bond.getAtom(1).getIndex());
            this.connections[bond.getAtom(1).getIndex()].add(bond.getAtom(0).getIndex());
        }
        for (int i=0; i<this.connections.length; i++) {
            if (this.connections[i].size() == 1) {
                tmpChainList.add(new ArrayList<>());
                tmpChainList.get(tmpChainList.size()-1).add(i);
            }
        }
        boolean tmpIsSearching = true;
        while (tmpIsSearching) {
            int tmpChainListIndex = 0;
            while (tmpChainListIndex < tmpChainList.size()) {
                List<Integer> tmpChainListAtIndex = tmpChainList.get(tmpChainListIndex);
                if (this.connections[tmpChainListAtIndex.get(tmpChainListAtIndex.size() - 1)].size() > 0 &&
                        tmpChainList.size() > 1) {
                    tmpChainListAtIndex.add(this.connections[tmpChainListAtIndex.get(tmpChainListAtIndex.size() - 1)].get(0));
                    this.connections[tmpChainListAtIndex.get(tmpChainListAtIndex.size() - 1)]
                            .remove(tmpChainListAtIndex.get(tmpChainListAtIndex.size() - 2));
                    if (this.connections[tmpChainListAtIndex.get(tmpChainListAtIndex.size() - 1)].size() > 1) {
                        this.branches.add(reverseList(tmpChainListAtIndex));
                        tmpChainList.remove(tmpChainListAtIndex);
                        tmpChainListIndex--;
                    }
                } else {
                    if (tmpChainList.size() == 2) {
                        this.branches.add(tmpChainList.get(0));
                        tmpChainList.remove(0);
                    }
                    this.branches.get(this.branches.size() - 1).remove(this.branches.get(this.branches.size() -1).size() - 1);
                    this.branches.get(this.branches.size() - 1).addAll(reverseList(tmpChainList.get(0)));
                    tmpChainList.remove(0);
                    tmpChainListIndex--;
                    tmpIsSearching = false;
                }
                tmpChainListIndex++;
            }
        }
    }

    private List<Integer> reverseList(List<Integer> aList) {
        List<Integer> revList = new ArrayList<>();
        for (int li : aList) {
            revList.add(0,li);
        }
        return revList;
    }

    /**
     * This method cuts linear molecule fragments from the this.branches list into chains of equal length which are
     * then stored in the this.fragmentsIndices list. Fragments that are too small or that
     */
    private void cutChains() {
        /*
        This method filters out all chain fragments that are smaller than minCut and stores them in the "rest" list.
        Likewise, it stores
        The first integer in each of these lists is the index of the atom that the chain is connected to.

        At first the method iterates through the chain list that consists of linear chain fragments.
         */
        this.rest = new ArrayList<>();
        this.fragmentsIndices = new ArrayList<>(this.molecule.getAtomCount()/this.minCut);
        List<Integer> tmpBranchingIndices = new ArrayList<>(this.branches.size()-1);
        for (List<Integer> tmpBranch : this.branches) {
            tmpBranchingIndices.add(tmpBranch.get(0));
        }
        int tmpBranchesIndex = 0;
        while (tmpBranchesIndex < this.branches.size()) {
            List<Integer> tmpBranchesItem = this.branches.get(tmpBranchesIndex);
            int tmpIndex = 1;
            int tmpBranchRest = 0;
            /*
            This while loop uses tmpIndex to iterate through the current tmpBranchesItem if it is not the last one in 
            the this.branches list. The last tmpBranchesItem is the longest chain in the molecule. All the others' first 
            integers are the indices of the atoms where the branches are connected to another molecule chain.
            If isPreservingTertiaryQuaternaryCarbons is true and/or a multiple bond would be cut, the tmpBranchRest 
            increases. Later the tmpBranchesRest will be used to remove the part of this branch, which needs be added
            back to the connected branch.
             */
            while (tmpIndex-tmpBranchRest < 2 && tmpIndex < tmpBranchesItem.size() && tmpBranchesIndex < this.branches.size()-1) {
                if (tmpIndex == 1 && isPreservingTertiaryQuaternaryCarbons) {
                    tmpBranchRest++;
                }
                else if (this.molecule.getBond(this.molecule.getAtom(tmpBranchesItem.get(tmpBranchRest)),
                        this.molecule.getAtom(tmpBranchesItem.get(tmpIndex))).getOrder() != IBond.Order.SINGLE) {
                    tmpBranchRest++;
                }
                tmpIndex++;
            }
            /*
            If the rest of the current branch (without the part that needs to be added back to the connected branch) is
            smaller than the minimum chain length (this.minCut) it means the branch is too small so that it will be
            added to the connected branch completely.
            If it is big enough though only the first part of the branch will be added back to the connected branch (to
            preserve tertiary and quaternary carbon atoms and to not split multiple bonds.
             */
            if (tmpBranchesItem.size() - tmpBranchRest <= this.minCut) {
                this.rest.add(tmpBranchesItem);
            } else {
                if (tmpBranchRest > 0 && tmpBranchesIndex < this.branches.size()-1) {
                    this.rest.add(tmpBranchesItem.subList(0, tmpBranchRest+1));
                }
                if (tmpBranchesIndex < this.branches.size()-1) {
                    tmpBranchesItem = tmpBranchesItem.subList(tmpBranchRest+1, tmpBranchesItem.size());
                }
                int tmpShift = 0;
                int tmpShift0 = 0;
                tmpIndex = 0;
                if (this.maxCut > 0) {
                    /*
                    If there is a set maximum chain length (this.maxCut > 0) the following while loop cuts the current
                    branch into chain fragments of equal length. In case a multiple bond would be cut, the current
                    fragment will be made smaller to shift the split to the previous bond. If the resulting fragment is
                    smaller than the minimum chain length the fragment will increase in size instead until the split
                    will be at a single bond.
                    All these fragments are then stored into the this.fragmentsIndices list.
                     */
                    while ((tmpIndex+1) * this.maxCut + tmpShift <= tmpBranchesItem.size()) {
                        tmpShift0 = tmpShift;
                        boolean tmpIsReversedShift = false;
                        if ((tmpIndex + 1) * this.maxCut + tmpShift + 1 < tmpBranchesItem.size()) {
                            while ((tmpIndex + 1) * this.maxCut + tmpShift + 1 < tmpBranchesItem.size() &&
                                    (this.molecule.getBond(this.molecule.getAtom(tmpBranchesItem.get((tmpIndex + 1) * this.maxCut + tmpShift)),
                                    this.molecule.getAtom(tmpBranchesItem.get((tmpIndex + 1) * this.maxCut + tmpShift + 1))).getOrder() != IBond.Order.SINGLE ||
                                    isPreservingTertiaryQuaternaryCarbons &&
                                    tmpBranchingIndices.contains(tmpBranchesItem.get((tmpIndex + 1) * this.maxCut + tmpShift)))) {
                                if (tmpShift0 - tmpShift < this.maxCut - this.minCut && !tmpIsReversedShift) {
                                    tmpShift--;
                                }
                                if (tmpShift0 - tmpShift >= this.maxCut - this.minCut && !tmpIsReversedShift) {
                                    tmpIsReversedShift = true;
                                    tmpShift = tmpShift0;
                                }
                                if (tmpIsReversedShift) {
                                    tmpShift++;
                                }
                            }
                        }
                        this.fragmentsIndices.add(tmpBranchesItem.subList(tmpIndex*this.maxCut+tmpShift0, (tmpIndex+1)*this.maxCut+tmpShift));
                        tmpIndex++;
                    }
                    /*
                    If the
                     */
                    if (tmpBranchesItem.size()%this.maxCut != 0) {
                        this.fragmentsIndices.add(tmpBranchesItem.subList(tmpIndex*this.maxCut+tmpShift0, tmpBranchesItem.size()));
                    }
                } else {
                    while ((tmpIndex + 1) * this.minCut + tmpShift <= tmpBranchesItem.size()) {
                        tmpShift0 = tmpShift;
                        if ((tmpIndex + 1) * this.minCut + tmpShift + 1 < tmpBranchesItem.size()) {
                            while ((tmpIndex + 1) * this.minCut + tmpShift + 1 < tmpBranchesItem.size() &&
                                    (this.molecule.getBond(this.molecule.getAtom(tmpBranchesItem.get((tmpIndex + 1) * this.minCut + tmpShift)),
                                    this.molecule.getAtom(tmpBranchesItem.get((tmpIndex + 1) * this.minCut + tmpShift + 1))).getOrder() != IBond.Order.SINGLE ||
                                    isPreservingTertiaryQuaternaryCarbons &&
                                    tmpBranchingIndices.contains(tmpBranchesItem.get((tmpIndex + 1) * this.minCut + tmpShift)))) {
                                tmpShift++;
                            }
                        }
                        this.fragmentsIndices.add(tmpBranchesItem.subList(tmpIndex * this.minCut + tmpShift0, (tmpIndex + 1) * this.minCut + tmpShift));
                        tmpIndex++;
                    }
                    if (tmpBranchesItem.size()%this.minCut != 0) {
                        this.rest.add(tmpBranchesItem.subList(tmpIndex*this.minCut+tmpShift0-1, tmpBranchesItem.size()));
                    }
                }
            }
            tmpBranchesIndex++;
        }
    }

    /**
     * During cutBranches all branches are separated from each other, but in order to preserve certain properties
     * through the fragmentation process the methods cutRings and cutChains produce rest fragments that need to be
     * added back to other branches, which is done by makeCorrections.
     * Every rest fragment has the index of its connecting atom in another branch stored. This
     */
    private void makeCorrections () {
        while (this.rest.size() > 0) {
            int tmpBranchesIndex = 0;
            while (tmpBranchesIndex < this.fragmentsIndices.size()) {
                List<Integer> tmpChainsAtIndex = this.fragmentsIndices.get(tmpBranchesIndex);
                int tmpRestIndex = 0;
                boolean tmpIsCombined = false;
                while (tmpRestIndex < this.rest.size() && !tmpIsCombined) {
                    List<Integer> restAtIndex = this.rest.get(tmpRestIndex);
                    if (tmpChainsAtIndex.contains(restAtIndex.get(0))) {
                        List<Integer> tmpCombinedFragments = new ArrayList<>(tmpChainsAtIndex);
                        tmpCombinedFragments.addAll(restAtIndex.subList(1, restAtIndex.size()));
                        this.fragmentsIndices.add(tmpCombinedFragments);
                        this.fragmentsIndices.remove(tmpChainsAtIndex);
                        this.rest.remove(restAtIndex);
                        tmpRestIndex++;
                        tmpIsCombined = true;
                    }
                    tmpRestIndex++;
                }
                if (!tmpIsCombined) {
                    tmpBranchesIndex++;
                }
            }
        }
    }
}