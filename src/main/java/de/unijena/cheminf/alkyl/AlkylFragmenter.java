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
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import java.util.ArrayList;
import java.util.List;

/**
 * The class AlkylFragmenter enables the user to dissect a hydrocarbon molecule into fragments of defined size and
 * properties. It separates each molecule branch (see cutBranches), and clips out ring linkers from ring clusters
 * (see cutRings). After splitting the molecule into its basic units, they are cut into fragments of the desired size,
 * which is between a minimum and a maximum size (see cutChains). Multiple bonds are always retained and the user can
 * choose whether also to preserve tertiary and quaternary carbon atoms, which means that after the fragmentation they
 * still have all their neighbouring atoms. These rules require parts of the branches to be added back to where they
 * were cut off (see makeCorrections).
 */
public class AlkylFragmenter {
    //<editor-fold desc="fragmentation settings">
    /**
     * IAtomContainer with all molecular information of the molecule that is to be fragmented. Variable is not changed
     * throughout the fragmentation process.
     */
    private IAtomContainer molecule;
    /**
     * User setting for the minimum fragment size.
     */
    private int minCut;
    /**
     * User setting for the maximum fragment size.
     */
    private int maxCut;
    /**
     * User setting for whether tertiary and quaternary carbon atoms are to be preserved during the fragmentation.
     */
    private boolean isPreservingTertiaryQuaternaryCarbons;
    //</editor-fold>
    //<editor-fold desc="private cache lists">
    /**
     * An array of Integer lists. Each list contains the indices of all neighbouring atoms of an atom. The array index
     * equals the atom's index in the IAtomContainer.
     */
    private List<Integer>[] connections;
    /**
     * A List of Integer Lists for each individual branch. It is filled by cutBranches and processed by cutChains.
     */
    private List<List<Integer>> branches;
    /**
     * List of Integer Lists for fragment remainders from the cutChains method. They will be added to their adjacent
     * branches in the makeCorrections method.
     */
    private List<List<Integer>> remainder;
    /**
     * List of Integer Lists for each resulting fragment. The cutRings and the cutChains method add to this list and the
     * makeCorrections method adds fragment remainders back to their adjacent fragments.
     */
    private List<List<Integer>> fragmentsIndices;
    /**
     * List of IAtomContainer objects with the same fragments as in the fragmentsIndices list. In addition to the atom
     * indices, they also contain bond information.
     */
    private List<IAtomContainer> fragmentsAtomContainer;
    //</editor-fold>
    //<editor-fold desc="public methods">
    /**
     * Public method for the user to commit a molecule to be fragmented.
     * @param aMolecule An IAtomContainer object containing the molecule to be fragmented.
     */
    public void setMolecule(IAtomContainer aMolecule) {
        this.molecule = aMolecule;
    }

    /**
     * Public method for the user to access the current molecule to be fragmented.
     * @return IAtomContainer object of the current molecule to be fragmented.
     */
    public IAtomContainer getMolecule() {
        return this.molecule;
    }
    /**
     * Public method for the user to set the fragmentation parameters and to start the fragmentation.
     * @param aMinCut Integer value for the minimal fragment size.
     * @param aMaxCut Integer value for the maximal fragment size.
     * @param aIsPreservingTertiaryQuaternaryCarbons Boolean for whether to preserve tertiary and quaternary carbons.
     * @throws CDKException Is triggered in genAtomContainer.
     */
    public void setFragmentationSettings(int aMinCut, int aMaxCut, boolean aIsPreservingTertiaryQuaternaryCarbons) throws CDKException {
        this.minCut = aMinCut;
        this.maxCut = aMaxCut;
        this.isPreservingTertiaryQuaternaryCarbons = aIsPreservingTertiaryQuaternaryCarbons;

        List<Integer> tmpAtomIndices = new ArrayList<>(this.molecule.getAtomCount());
        for (int i=0; i<this.molecule.getAtomCount(); i++) {
            tmpAtomIndices.add(i);
        }
        cutBranches(tmpAtomIndices);
        cutChains();
        cutRings();
        makeCorrections();
        genAtomContainer(this.fragmentsIndices);
    }

    /**
     * Public method for the user to access the resulting fragment molecules as IAtomContainer objects.
     * @return An ArraryList of IAtomContainers containing the resulting fragment molecules.
     */
    public List<IAtomContainer> getIAtomContainer () {
        return this.fragmentsAtomContainer;
    }
    //</editor-fold>

    //<editor-fold desc="test methods">
    /**
     * Method to examine the results from the cutBranches method.
     * @return An ArrayList of Integer ArrayList objects containing the atom indices of each molecular branch.
     */
    protected List<List<Integer>> getBranches () {
        List<Integer> tmpAtomIndices = new ArrayList<>(this.molecule.getAtomCount());
        for (int i=0; i<this.molecule.getAtomCount(); i++) {
            tmpAtomIndices.add(i);
        }
        cutBranches(tmpAtomIndices);
        return this.branches;
    }
    //</editor-fold>
    //<editor-fold desc="private fragmentation methods">

    /**
     * This method dissects a given part of the molecule into its individual branches.
     * @param anBranchedMoleculeFragment An Integer ArrayList of atom indices of molecule fragment.
     */
    private void cutBranches (List<Integer> anBranchedMoleculeFragment) {
        /*
        At the beginning an array of lists is created. Each array item represents an atom with its neighbouring atoms 
        stored in a list. 
         */
        List<Integer>[] tmpConnections = new ArrayList[this.molecule.getAtomCount()];
        for (int i = 0; i < tmpConnections.length; i++) {
            tmpConnections[i] = new ArrayList<>();
        }
        for (IBond bond : this.molecule.bonds()) {
            int atom0 = bond.getAtom(0).getIndex();
            int atom1 = bond.getAtom(1).getIndex();
            if (anBranchedMoleculeFragment.contains(atom0) && anBranchedMoleculeFragment.contains(atom1)) {
                tmpConnections[atom0].add(atom1);
                tmpConnections[atom1].add(atom0);
            }
        }
        List<List<Integer>> tmpChainList = new ArrayList<>();
        
        /*
        The indices of all primary (terminal) carbon atoms become each the start of a new chain (beginning of new list 
        in tmpChainList).  
         */
        for (int i=0; i<tmpConnections.length; i++) {
            if (tmpConnections[i].size() == 1) {
                tmpChainList.add(new ArrayList<>());
                tmpChainList.get(tmpChainList.size()-1).add(i);
            }
        }
        /*
        Next, all chains grow atom by atom at the same time; meaning, the first chain in tmpChainList adds its 
        neighbouring atom, then the second chain gets its adjacent atom added and so on. After all chains have added 
        their adjacent atoms they move on to the next. Every time one chain encounters a branching it is removed from
        the tmpChainList and added to the this.branches list. At the end, there are only two chains left which will meet
        at the centre of the molecule and then joint together to become the main branch. At this point, the variable 
        tmpIsSearing becomes false and the loop finishes.
         */
        if (this.branches == null) {
            this.branches = new ArrayList<>();
            this.connections = tmpConnections;
            this.remainder = new ArrayList<>();
            if (this.minCut > 0) {
                this.fragmentsIndices = new ArrayList<>(this.molecule.getAtomCount() / this.minCut);
            } else {
                this.fragmentsIndices = new ArrayList<>(this.branches.size());
            }
        }
        boolean tmpIsSearching = true;
        while (tmpIsSearching) {
            int tmpChainListIndex = 0;
            while (tmpChainListIndex < tmpChainList.size()) {
                List<Integer> tmpCurrentChain = tmpChainList.get(tmpChainListIndex);
                if (tmpConnections[tmpCurrentChain.get(tmpCurrentChain.size() - 1)].size() > 0) {
                    /*
                    If tmpConnections[tmpCurrentChain.get(tmpCurrentChain.size() - 1)].size() is zero it means that
                    the last atom in the current chain does not have neighbouring atoms. This happens when the two last
                    chains have reached each other.
                    Before that happens, the neighbouring atom is added to the current chain and the bond information
                    erased from the tmpConnections array. If the newly added atom has more than one other neighbour
                    (at a branching) the whole current chain is added to the this.branches list and removed from the
                    tmpChainList. This way each subbranch also contains the information, from which atom of another
                    branch it was cut off.
                     */
                    tmpCurrentChain.add(tmpConnections[tmpCurrentChain.get(tmpCurrentChain.size() - 1)].get(0));
                    tmpConnections[tmpCurrentChain.get(tmpCurrentChain.size() - 1)]
                            .remove(tmpCurrentChain.get(tmpCurrentChain.size() - 2));
                    tmpConnections[tmpCurrentChain.get(tmpCurrentChain.size() - 2)]
                            .remove(tmpCurrentChain.get(tmpCurrentChain.size() - 1));
                    if (tmpConnections[tmpCurrentChain.get(tmpCurrentChain.size() - 1)].size() > 1) {
                        this.branches.add(reverseList(tmpCurrentChain));
                        tmpChainList.remove(tmpCurrentChain);
                        tmpChainListIndex--;
                    }
                } else {
                    /*
                    When the two last chains meet each other, the first chain is added to the this.branches list and the
                    second is then merged together with the just added chain.
                     */
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
            /*
            When all branches are cut off, but if the last two chains do not meet each other it means that there must be
            a ring (cluster) in between, which will be further broken down by the cutRings method. In this case an empty
            list will be added at the end of the this.branches list to indicate that there is no continuous main branch.
             */
            if (tmpIsSearching && tmpChainList.size() == 0) {
                this.branches.add(new ArrayList<>());
                tmpIsSearching = false;
            }

        }
    }

    /**
     * The reverseList method reverses ArrayLists of Integers.
     * @param aList an Integer ArrayList
     * @return reversed Integer ArrayList
     */
    private List<Integer> reverseList(List<Integer> aList) {
        List<Integer> revList = new ArrayList<>();
        for (int li : aList) {
            revList.add(0,li);
        }
        return revList;
    }

    /**
     * After the cutBranches method extracted only the non-cyclic parts of the molecule, the remaining rings and ring
     * linkers are separated from each other by the cutRings method.
     */
    private void cutRings () {
        List<List<Integer>> tmpFragment = new ArrayList<>();
        List<Integer> tmpCurrentChain = new ArrayList<>();
        List<Integer> tmpBranchStarter = new ArrayList<>();
        List<Integer> tmpFragmentStarter = new ArrayList<>();
        boolean tmpIsRingNotRingLinker;
        boolean tmpHasSameKindNeighbour;
        Integer tmpIndexCurrentAtom;
        Cycles.markRingAtomsAndBonds(molecule);
        /*
        The cutBranches method extracted all non-cyclic and non-ring-linker atoms from the this.connections array. At
        the beginning of the cutRings method the first ring or ring linker atom from the molecule is being identified by
        checking which of the atoms still has an entry about neighbouring atoms. This atom is saved into the
        tmpFragmentStarters list and is going to be the starting atom for the iteration.
         */
        int tmpConnectionsIndex = 0;
        while (tmpFragmentStarter.size() == 0 && tmpConnectionsIndex < this.connections.length) {
            if (this.connections[tmpConnectionsIndex].size() != 0) {
                tmpFragmentStarter.add(tmpConnectionsIndex);
            }
            tmpConnectionsIndex++;
        }
        while (tmpFragmentStarter.size() != 0) {
            /*
            Each atom in the tmpFragmentStarter is the beginning of either a ring or a ring linker fragment.
            They become the start of new chains to follow (tmpBranchStarter).
             */
            tmpIndexCurrentAtom = tmpFragmentStarter.get(0);
            tmpFragmentStarter.remove(0);
            tmpIsRingNotRingLinker = this.molecule.getAtom(tmpIndexCurrentAtom).isInRing();
            tmpBranchStarter.add(tmpIndexCurrentAtom);
            while (tmpBranchStarter.size() != 0) {
                /*
                This loop looks for atoms of the same type (ring or ring linker atom) and adds them either to the
                tmpCurrentChain or the tmpBranchStarter if there are more than one atom of the same kind at a branching,
                until the fragment is finished. All other atoms are added to the tmpFragmentStarter list to start new
                fragments.

                But there is a difference of procedure when tertiary and quaternary carbon atoms need to be preserved.
                Then the
                */
                tmpIndexCurrentAtom = tmpBranchStarter.get(0);
                if (!tmpCurrentChain.contains(tmpIndexCurrentAtom)) {
                    tmpCurrentChain.add(tmpIndexCurrentAtom);
                }
                tmpBranchStarter.remove(0);

                if (this.isPreservingTertiaryQuaternaryCarbons) {
                    /*
                    When starting with a new tmpCurrentChain, there might be more than one neighbouring atom to
                    tmpIndexCurrentAtom. The first becomes tmpIndexNeighbouringAtom and all the others become branch
                    starters.
                     */
                    while (this.connections[tmpIndexCurrentAtom].size() > 0) {
                        /*
                        This loop continues until it reaches a previous atom. Because all other neighbouring atoms of
                        tmpIndexCurrentAtom have been added to tmpBranchStarter or tmpFragmentStarter only one
                        neighbouring atom remains.
                         */
                        int tmpIndexNeighbouringAtom = this.connections[tmpIndexCurrentAtom].get(0);
                        this.connections[tmpIndexNeighbouringAtom].remove((Integer) tmpIndexCurrentAtom);
                        for (int tmpIndex=1; tmpIndex<this.connections[tmpIndexCurrentAtom].size(); tmpIndex++) {
                            tmpBranchStarter.add(this.connections[tmpIndexCurrentAtom].get(tmpIndex));
                            this.connections[tmpBranchStarter.get(tmpBranchStarter.size()-1)].remove((Integer) tmpIndexCurrentAtom);
                        }
                        int tmpIndexNextNeighbouringAtom = 0;
                        for (int tmpIndex=this.connections[tmpIndexNeighbouringAtom].size()-1; tmpIndex>=0; tmpIndex--) {
                            /*
                            Loop to access all tmpIndexNextNeighbouringAtom (neighbouring atoms to the neighbouring atom
                            of tmpIndexCurrentAtom).
                             */
                            tmpIndexNextNeighbouringAtom = this.connections[tmpIndexNeighbouringAtom].get(tmpIndex);
                            this.connections[tmpIndexNextNeighbouringAtom].remove((Integer) tmpIndexNeighbouringAtom);
                            if (!tmpIsRingNotRingLinker && (this.molecule.getAtom(tmpIndexNeighbouringAtom).isInRing() ||
                                    this.molecule.getAtom(tmpIndexNextNeighbouringAtom).isInRing() &&
                                            (this.molecule.getAtom(tmpIndexCurrentAtom).isInRing() ||
                                                    this.molecule.getAtom(tmpIndexNeighbouringAtom).isInRing()))) {
                                /*
                                If tmpIndexNeighbouringAtom or tmpIndexNextNeighbouringAtom is a ring atom while the
                                current atom type is ring linkers, the current type is changed to ring atoms.
                                 */
                                tmpIsRingNotRingLinker = true;
                            }
                            if (tmpIndex == 0 && this.connections[tmpIndexNeighbouringAtom].size() == 1 &&
                                    !this.molecule.getAtom(tmpIndexNeighbouringAtom).isInRing() &&
                                    !this.molecule.getAtom(tmpIndexNextNeighbouringAtom).isInRing()) {
                                    /*
                                    If both, the tmpIndexNeighbouringAtom and the tmpIndexNextNeighbouringAtom are ring
                                    linkers and not branched, the current branch ends with the tmpIndexNeighbouringAtom
                                    and the tmpIndexNextNeighbouringAtom is added to tmpFragmentStarter
                                     */
                                this.connections[tmpIndexCurrentAtom] = new ArrayList<>();
                                tmpFragmentStarter.add(tmpIndexNextNeighbouringAtom);
                                this.connections[tmpIndexNextNeighbouringAtom].remove((Integer) tmpIndexNeighbouringAtom);
                                tmpCurrentChain.add(tmpIndexNeighbouringAtom);
                                if (!tmpIsRingNotRingLinker) {
                                    tmpFragment.add(tmpCurrentChain);
                                    tmpCurrentChain = new ArrayList<>();
                                }
                            } else if (tmpIndex == 0 && this.connections[tmpIndexNeighbouringAtom].size() == 1 &&
                                    !this.molecule.getAtom(tmpIndexCurrentAtom).isInRing() &&
                                    !this.molecule.getAtom(tmpIndexNeighbouringAtom).isInRing()) {
                                this.connections[tmpIndexCurrentAtom] = new ArrayList<>();
                                tmpFragmentStarter.add(tmpIndexNeighbouringAtom);
                                if (!tmpIsRingNotRingLinker) {
                                    tmpFragment.add(tmpCurrentChain);
                                    tmpCurrentChain = new ArrayList<>();
                                }
                            } else if (tmpIndex > 0) {
                                tmpBranchStarter.add(tmpIndexNextNeighbouringAtom);
                            } else {
                                this.connections[tmpIndexNeighbouringAtom].remove((Integer) tmpIndexCurrentAtom);
                                tmpCurrentChain.add(tmpIndexNeighbouringAtom);
                                this.connections[tmpIndexCurrentAtom] = new ArrayList<>();
                                tmpIndexCurrentAtom = tmpIndexNeighbouringAtom;
                                tmpIndexNeighbouringAtom = tmpIndexNextNeighbouringAtom;
                            }
                        }
                        if (this.connections[tmpIndexNeighbouringAtom].size() == 0) {
                            tmpCurrentChain.add(tmpIndexNeighbouringAtom);
                            this.connections[tmpIndexCurrentAtom] = new ArrayList<>();
                        }
                    }
                } else {
                    /*
                    The first neighbouring atom that is of the same kind as the current atom (both ring or both
                    ring linker atoms) is added to the tmpCurrentChain list. All other neighbouring atoms of the same
                    kind are added to tmpBranchStarter. Those neighbouring atoms that are of the other kind are added to
                    the tmpFragmentStarter list.
                     */

                    tmpHasSameKindNeighbour = false;
                    while (this.connections[tmpIndexCurrentAtom].size() > 0) {
                        /*
                        This loop continues until the algorithm arrives at an earlier atom again. Then it starts again
                        with a new branch or a fragment starter (outer loops).
                         */
                        for (int tmpNeighbouringAtom : this.connections[tmpIndexCurrentAtom]) {
                            /*
                            In order to sort the atoms into rings and ring linkers, the method iterates through each of
                            the neighbours of every atom.
                             */
                            if ((this.molecule.getAtom(tmpNeighbouringAtom).isInRing() == tmpIsRingNotRingLinker ||
                                    !this.molecule.getAtom(tmpNeighbouringAtom).isInRing() == !tmpIsRingNotRingLinker)
                                    && !tmpHasSameKindNeighbour) {
                                /*
                                If the neighbour is the first atom of the same type (ring or ring linker), it is added
                                to the tmpCurrentChain list.
                                 */
                                tmpCurrentChain.add(tmpNeighbouringAtom);
                                tmpHasSameKindNeighbour = true;
                            } else if (this.molecule.getAtom(tmpNeighbouringAtom).isInRing() == tmpIsRingNotRingLinker
                                    && tmpHasSameKindNeighbour) {
                                /*
                                All other neighbours of the same type are added to the tmpBranchStarter list to be added
                                to the current fragment later on.
                                 */
                                tmpBranchStarter.add(tmpNeighbouringAtom);
                            } else if (this.molecule.getAtom(tmpNeighbouringAtom).isInRing() != tmpIsRingNotRingLinker) {
                                /*
                                All neighbouring atoms of another type are added to the tmpFragmentStarter list to
                                become their own fragments later on.
                                 */
                                tmpFragmentStarter.add(tmpNeighbouringAtom);
                            }
                            if (tmpIsRingNotRingLinker || !tmpHasSameKindNeighbour) {
                                /*
                                To make sure that the tmpCurrentChain does not go back to previous atoms, the
                                tmpIndexCurrentAtom is removed from the tmpNeighbouringAtom.
                                 */
                                this.connections[tmpNeighbouringAtom].remove((Integer) tmpIndexCurrentAtom);
                            }
                        }
                        /*
                        After all neighbouring atoms were placed in their respective list, the neighbours' list of the
                        tmpIndexCurrentAtom is cleared.
                         */
                        this.connections[tmpIndexCurrentAtom] = new ArrayList<>();
                        if (tmpHasSameKindNeighbour) {
                            /*
                            If the tmpIndexCurrentAtom has at least one neighbour of the same type, the next
                            tmpIndexCurrentAtom is the neighbouring atom, which was recently added to tmpCurrentChain.
                             */
                            tmpIndexCurrentAtom = tmpCurrentChain.get(tmpCurrentChain.size()-1);
                            tmpHasSameKindNeighbour = false;
                        } else if (!tmpHasSameKindNeighbour && tmpBranchStarter.size() > 0) {
                            /*
                            If all neighbouring atoms were of the other type, and if the tmpBranchStarter list is not
                            empty, the first atom in the list becomes the tmpIndexCurrentAtom.
                             */
                            tmpIndexCurrentAtom = tmpBranchStarter.get(0);
                            if (tmpCurrentChain.contains(tmpIndexCurrentAtom)) {
                                tmpBranchStarter.remove((Integer) tmpIndexCurrentAtom);
                            }
                        } else if (!tmpHasSameKindNeighbour && tmpFragmentStarter.size() > 0) {
                            /*
                            If all neighbouring atoms were of the other type, and if the tmpBranchStarter list is empty,
                            but not tmpFragmentStarter, the first atom in the list becomes the tmpIndexCurrentAtom.
                             */
                            tmpIndexCurrentAtom = tmpFragmentStarter.get(0);
                            if (tmpCurrentChain.contains(tmpIndexCurrentAtom)) {
                                tmpFragmentStarter.remove((Integer) tmpIndexCurrentAtom);
                            }
                        }
                    }
                }
            }
            if (tmpIsRingNotRingLinker) {
                /*
                If the current atom type is ring atom the fragment will not be further broken down and it is added to
                the this.fragmentsIndices list.
                 */
                this.fragmentsIndices.add(tmpCurrentChain);
            } else {
                /*
                Ring Linker fragments are first separated into individual branches in cutBranches and cut in the desired
                size in cutChains.
                 */
                cutBranches(tmpCurrentChain);
                cutChains();
            }
            tmpCurrentChain = new ArrayList<>();
        }
        /*
        The last part of this method is only used when tmpIsPreservingTertiaryQuaternaryCarbon is enabled. It iterates
        through the tmpFragment list and looks for those fragments that belong together by checking if the neighbouring
        atom of the last atom in the fragment list is also contained in another fragment list.
         */
        while (tmpFragment.size() > 0) {
            List<Integer> tmpMergedFragments = tmpFragment.get(0);
            boolean tmpListHasChanged = true;
            while (tmpListHasChanged) {
                tmpListHasChanged = false;
                int tmpFragmentIndex = 1;
                while (tmpFragmentIndex < tmpFragment.size()){
                    List<Integer> tmpCurrentFragment = tmpFragment.get(tmpFragmentIndex);
                    if (this.connections[tmpCurrentFragment.get(tmpCurrentFragment.size()-1)].size() > 0 &&
                            tmpMergedFragments.contains(this.connections[tmpCurrentFragment.get(tmpCurrentFragment.size()-1)].get(0)) ||
                            this.connections[tmpFragment.get(0).get(tmpFragment.get(0).size()-1)].size() > 0 &&
                                    tmpCurrentFragment.contains(this.connections[tmpFragment.get(0).get(tmpFragment.get(0).size()-1)].get(0))) {
                        tmpMergedFragments.addAll(tmpCurrentFragment);
                        tmpFragment.remove(tmpCurrentFragment);
                        tmpFragmentIndex--;
                        tmpListHasChanged = true;
                    }
                    tmpFragmentIndex++;
                }
            }
            cutBranches(tmpMergedFragments);
            cutChains();
            tmpFragment.remove(0);
        }
    }

    /**
     * This method cuts linear molecule fragments from the this.branches list into chains of equal length which are
     * then stored in the this.fragmentsIndices list. It does not cut multiple bonds and 
     */
    private void cutChains() {
        /*
        The tmpBranchingIndices list is used for the tmpIsPreservingTertiaryQuaternaryCarbon option. It contains the
        indices of all branching atoms in a chain.
         */
        List<Integer> tmpBranchingIndices = new ArrayList<>(this.branches.size()-1);
        if (this.branches.get(0).size() > 0) {
            for (List<Integer> tmpBranch : this.branches) {
                if (!tmpBranchingIndices.contains(tmpBranch.get(0))) {
                    tmpBranchingIndices.add(tmpBranch.get(0));
                }
            }
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
            while (tmpIndex - tmpBranchRest < 2 && tmpIndex < tmpBranchesItem.size()
                    && tmpBranchesIndex < this.branches.size()-1) {
                if (tmpIndex == 1 && isPreservingTertiaryQuaternaryCarbons ||
                        this.molecule.getBond(this.molecule.getAtom(tmpBranchesItem.get(tmpBranchRest)),
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
            if (tmpBranchesItem.size() - tmpBranchRest <= this.minCut && tmpBranchesItem.size() - tmpBranchRest > 0) {
                this.remainder.add(tmpBranchesItem);
            } else {
                if (tmpBranchRest > 0 && tmpBranchesIndex < this.branches.size()-1) {
                    this.remainder.add(tmpBranchesItem.subList(0, tmpBranchRest+1));
                }
                if (tmpBranchesIndex < this.branches.size()-1) {
                    tmpBranchesItem = tmpBranchesItem.subList(tmpBranchRest+1, tmpBranchesItem.size());
                }
                if (this.minCut == 0 && this.maxCut == 0 && tmpBranchesItem.size() > 0) {
                    this.fragmentsIndices.add(tmpBranchesItem);
                }
                /*
                tmpIndexCutPosition stands for the index of the bond within the tmpBranchesItem. E.g. if it is 3, the
                third bond will be broken. tmpIndexCutPosition is the index of the next cut position. Because multiple
                bonds must not be split, there are two variables: tmpShift and tmpShift0. They shift or move the cut
                for- or backward along the chain.
                 */
                int tmpIndexCutPosition;
                int tmpIndexNextCutPosition;
                int tmpShift = 0;
                int tmpShift0 = 0;
                if (this.maxCut > 0) {
                    /*
                    If there is a set maximum chain length (this.maxCut > 0) the following while loop cuts the current
                    branch into chain fragments of equal length. In case a multiple bond would be cut, the current
                    fragment will be made smaller to shift the split to the previous bond. If the resulting fragment is
                    smaller than the minimum chain length the fragment will increase in size instead until the split
                    will be at a single bond.
                    All these fragments are then stored into the this.fragmentsIndices list.
                     */
                    tmpIndexCutPosition = tmpBranchesItem.size();
                    tmpIndexNextCutPosition = tmpIndexCutPosition-this.maxCut;
                    while (tmpIndexNextCutPosition + tmpShift >= 0) {
                        tmpShift0 = tmpShift;
                        boolean tmpIsReversedShift = false;
                        while (tmpIndexNextCutPosition + tmpShift - 1 > 0 &&
                                (this.molecule.getBond(this.molecule.getAtom(tmpBranchesItem.get(tmpIndexNextCutPosition+tmpShift)),
                                this.molecule.getAtom(tmpBranchesItem.get(tmpIndexNextCutPosition+tmpShift-1))).getOrder() != IBond.Order.SINGLE ||
                                (tmpBranchingIndices.contains(tmpBranchesItem.get(tmpIndexNextCutPosition+tmpShift)) ||
                                        tmpBranchingIndices.contains(tmpBranchesItem.get(tmpIndexNextCutPosition+tmpShift-1)))
                                        && isPreservingTertiaryQuaternaryCarbons)) {
                            /*
                            The conditions are: Continue as long as 1. tmpIndexNextCutPosition is within the chain
                            length range, 2. there is a multiple bond between the current and the next atom, or 3. the
                            current atom is next to a branching when isPreservingTertiaryQuaternaryCarbons is enabled.
                             */
                            if (tmpShift - tmpShift0 < this.maxCut - this.minCut && !tmpIsReversedShift) {
                                tmpShift++;
                            }
                            if (tmpShift - tmpShift0 >= this.maxCut - this.minCut && !tmpIsReversedShift) {
                                tmpIsReversedShift = true;
                                tmpShift = tmpShift0;
                            }
                            if (tmpIsReversedShift) {
                                tmpShift--;
                            }
                        }
                        /*
                        After having determined the right tmpShift value, the embedded fragment piece is added to the
                        this.fragmentsIndices list.
                         */
                        this.fragmentsIndices.add(tmpBranchesItem.subList(tmpIndexNextCutPosition+tmpShift, tmpIndexCutPosition+tmpShift0));
                        tmpIndexCutPosition -= this.maxCut;
                        tmpIndexNextCutPosition -= this.maxCut;
                    }
                    /*
                    If the chain cannot be divided evenly without a remainder, what is left is either also added to
                    this.fragmentsIndices or to this.remainders if it is smaller than this.minCut.
                     */
                    if (tmpIndexCutPosition+tmpShift % this.maxCut != 0) {
                        if (tmpIndexCutPosition+tmpShift >= this.minCut) {
                            this.fragmentsIndices.add(tmpBranchesItem.subList(0, tmpIndexCutPosition + tmpShift));
                        } else {
                            this.remainder.add(reverseList(tmpBranchesItem.subList(0, tmpIndexCutPosition + tmpShift + 1)));
                        }
                    }
                } else if (this.maxCut == 0 && this.minCut > 0) {
                    /*
                    If there is only a this.minCut, the whole tmpShift finding process becomes simpler. It only shifts
                    into one direction. But otherwise this part of the method is the same as for the this.maxCut part.
                     */
                    tmpIndexCutPosition = tmpBranchesItem.size();
                    tmpIndexNextCutPosition = tmpIndexCutPosition-this.minCut;
                    while (tmpIndexNextCutPosition + tmpShift >= 0) {
                        tmpShift0 = tmpShift;
                        while (tmpIndexNextCutPosition + tmpShift - 1 > 0 &&
                                (this.molecule.getBond(this.molecule.getAtom(tmpBranchesItem.get(tmpIndexNextCutPosition + tmpShift)),
                                this.molecule.getAtom(tmpBranchesItem.get(tmpIndexNextCutPosition+tmpShift-1))).getOrder() != IBond.Order.SINGLE ||
                                (tmpBranchingIndices.contains(tmpBranchesItem.get(tmpIndexNextCutPosition))) && isPreservingTertiaryQuaternaryCarbons)) {
                            tmpShift--;
                        }
                        this.fragmentsIndices.add(tmpBranchesItem.subList(tmpIndexNextCutPosition+tmpShift, tmpIndexCutPosition+tmpShift0));
                        tmpIndexCutPosition -= this.minCut;
                        tmpIndexNextCutPosition -= this.minCut;
                    }
                    if (tmpIndexCutPosition+tmpShift % this.minCut != 0) {
                        if (tmpIndexCutPosition+tmpShift >= this.minCut) {
                            this.fragmentsIndices.add(tmpBranchesItem.subList(0, tmpIndexCutPosition + tmpShift));
                        } else {
                            this.remainder.add(reverseList(tmpBranchesItem.subList(0, tmpIndexCutPosition + tmpShift + 1)));
                        }
                    }
                }
            }
            this.branches.remove(tmpBranchesIndex);
        }
    }

    /**
     * During cutBranches all branches are separated from each other, but in order to preserve certain properties
     * through the fragmentation process the methods cutRings and cutChains produce rest fragments that need to be
     * added back to other branches, which is done by makeCorrections.
     * Every remainder fragment also contains the index of its connecting atom in the adjacent branch.
     */
    private void makeCorrections () {
        while (this.remainder.size() > 0) {
            int tmpBranchesIndex = 0;
            while (tmpBranchesIndex < this.fragmentsIndices.size()) {
                List<Integer> tmpChainsAtIndex = this.fragmentsIndices.get(tmpBranchesIndex);
                int tmpRestIndex = 0;
                boolean tmpIsCombined = false;
                while (tmpRestIndex < this.remainder.size() && !tmpIsCombined) {
                    List<Integer> restAtIndex = this.remainder.get(tmpRestIndex);
                    if (tmpChainsAtIndex.contains(restAtIndex.get(0))) {
                        List<Integer> tmpCombinedFragments = new ArrayList<>(tmpChainsAtIndex);
                        tmpCombinedFragments.addAll(restAtIndex.subList(1, restAtIndex.size()));
                        this.fragmentsIndices.add(tmpCombinedFragments);
                        this.fragmentsIndices.remove(tmpChainsAtIndex);
                        this.remainder.remove(restAtIndex);
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
    /**
     * This method converts lists of atom indices into IAtomContainer objects. Because during the fragmentation
     * algorithm only atom indices are used, the lists of fragments consisting of indices need to be converted back into
     * IAtomContainer objects. Because bonds are broken during the fragmentation, the fragments also need to be
     * saturated with hydrogen atoms.
     * @param anIndicesList An ArrayList with ArrayList objects containing the atom indices of fragment molecules of the
     *                      molecule to be fragmented.
     * @throws CDKException Is triggered when saturating fragment molecules with hydrogen where bonds were split.
     */

    private void genAtomContainer(List<List<Integer>> anIndicesList) throws CDKException {
        this.fragmentsAtomContainer = new ArrayList<>(anIndicesList.size());
        for (List<Integer> tmpListItem : anIndicesList) {
            IAtomContainer tmpMoleculeFragment = new AtomContainer();
            for (IBond tmpBond : molecule.bonds()) {
                if (tmpListItem.contains(tmpBond.getAtom(0).getIndex()) &&
                        tmpListItem.contains(tmpBond.getAtom(1).getIndex())) {
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
    //</editor-fold>
}