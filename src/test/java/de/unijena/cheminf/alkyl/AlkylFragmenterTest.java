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

package de.unijena.cheminf.alkyl;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 */
public class AlkylFragmenterTest extends AlkylFragmenter {
    //<editor-fold desc="Definition and Declaration of Private Objects">
    /**
     *
     */
    private SmilesParser sp;
    private AlkylFragmenter fragmenter;
    private Logger logger;
    public AlkylFragmenterTest () {
        this.sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        this.fragmenter = new AlkylFragmenter();
        this.logger = Logger.getLogger(AlkylFragmenterTest.class.getName());
    }
    //</editor-fold>
    //<editor-fold desc="">
    /**
     * Converts an ArrayList of IAtomContainer objects into an ArrayList of SMILES String objects.
     * @param aFragmentsList An ArrayList of IAtomContainers of molecule fragments.
     * @return An ArrayList of SMILES String objects of molecule fragments.
     * @throws CDKException Is triggered when IAtomContainer of a fragment is not convertible into a String.
     */
    private List<String> getSmiles(List<IAtomContainer> aFragmentsList) throws CDKException {
        SmilesGenerator tmpSmilesGenerator = SmilesGenerator.generic();
        List<String> tmpSmilesList = new ArrayList<>(aFragmentsList.size());
        for (IAtomContainer tmpFragment : aFragmentsList) {
            tmpSmilesList.add(tmpSmilesGenerator.create(tmpFragment));
        }
        return tmpSmilesList;
    }
    private String expectedActualLists (String[] aSmilesListExpected, List<String> aSmilesListActual) {
        String tmpOutput = "\nExpected: ";
        for (int i=0; i<aSmilesListExpected.length; i++) {
            tmpOutput += aSmilesListExpected[i];
            if (i < aSmilesListExpected.length-1) {
                tmpOutput += ", ";
            }
        }
        tmpOutput += "\nActual:   ";
        for (int i=0; i<aSmilesListActual.size(); i++) {
            tmpOutput += aSmilesListActual.get(i);
            if (i < aSmilesListActual.size()-1) {
                tmpOutput += ", ";
            }
        }
        return tmpOutput;
    }

    private String expectedActualLists (List<Integer>[] anIntegerListExpected, List<List<Integer>> anIntegerListActual) {
        String tmpOutput = "\nExpected: ";
        for (int i=0; i<anIntegerListExpected.length; i++) {
            tmpOutput += anIntegerListExpected[i];
            if (i < anIntegerListExpected.length-1) {
                tmpOutput += ", ";
            }
        }
        tmpOutput += "\nActual:   ";
        for (int i=0; i<anIntegerListActual.size(); i++) {
            tmpOutput += anIntegerListActual.get(i);
            if (i < anIntegerListActual.size()-1) {
                tmpOutput += ", ";
            }
        }
        return tmpOutput;
    }

    private boolean compareLists (String[] aSmilesListExpected, List<String> aSmilesListActual) {
        List<String> tmpSmilesList1 = new ArrayList<>(aSmilesListActual);
        if (aSmilesListExpected.length != aSmilesListActual.size()) {
            String tmpOutput = expectedActualLists(aSmilesListExpected, aSmilesListActual);
            tmpOutput += "\nThe expected and actual SMILES lists do not have the same size";
            logger.log(Level.INFO, tmpOutput);
            return false;
        }
        for (String tmpSmiles : aSmilesListExpected) {
            if (tmpSmilesList1.contains(tmpSmiles)) {
                tmpSmilesList1.remove(tmpSmiles);
            } else {
                String tmpOutput = expectedActualLists(aSmilesListExpected, aSmilesListActual);
                tmpOutput += "\nThe expected SMILES list does not contain " + tmpSmiles;
                logger.log(Level.INFO, tmpOutput);
                return false;
            }
        }
        return true;
    }

    private boolean compareLists (List<Integer>[] anIntegerListExpected, List<List<Integer>> anIntegerListActual) {
        List<List<Integer>> tmpSmilesList1 = anIntegerListActual;
        if (anIntegerListExpected.length != anIntegerListActual.size()) {
            String tmpOutput = expectedActualLists(anIntegerListExpected, anIntegerListActual);
            tmpOutput += "\nThe expected and actual list of atom indices do not have the same size";
            logger.log(Level.INFO, tmpOutput);
            return false;
        }
        for (List<Integer> tmpSmiles : anIntegerListExpected) {
            if (tmpSmilesList1.contains(tmpSmiles)) {
                tmpSmilesList1.remove(tmpSmiles);
            } else {
                String tmpOutput = expectedActualLists(anIntegerListExpected, anIntegerListActual);
                tmpOutput += "\nThe expected list of atom indices does not contain " + tmpSmiles;
                logger.log(Level.INFO, tmpOutput);
                return false;
            }
        }
        return true;
    }

    private void testNumberOfAtoms () {
        int tmpNumberOfAtomsInFragmentList = 0;
        for (IAtomContainer tmpFragment : this.fragmenter.getIAtomContainer()) {
            tmpNumberOfAtomsInFragmentList += tmpFragment.getAtomCount();
        }
        Assert.assertEquals(this.fragmenter.getMolecule().getAtomCount(), tmpNumberOfAtomsInFragmentList);
    }
    //</editor-fold>
    //<editor-fold desc="Test Fragmentation Methods">
    //<editor-fold desc="cutBranches">
    /**
     * Test method to examine the cutBranches method for non-cyclic molecules.
     * @throws InvalidSmilesException Is triggered when SMILES code does not translate into a molecule.
     */
    @Test
    public void testCutBranchesNonCyclic () throws InvalidSmilesException {
        IAtomContainer mol1 = this.sp.parseSmiles("CCCCC(CCC)C(C)(CC)C(C)CC(C)(C)C");
        this.fragmenter.setMolecule(mol1);
        Assert.assertTrue(compareLists(new List[]{List.of(8,9),List.of(12,13),List.of(15,16),List.of(15,17),List.of(8,10,11),List.of(4,5,6,7),List.of(0,1,2,3,4,8,12,14,15,18)},
                this.fragmenter.getBranches()));
    }
    /**
     * Test method to examine the cutBranches method for cyclic molecules.
     * @throws InvalidSmilesException Is triggered when SMILES code does not translate into a molecule.
     */
    @Test
    public void testCutBranchesCyclic () throws InvalidSmilesException {
        IAtomContainer mol1 = this.sp.parseSmiles("CCC(C)C1=CC(C)(CC)C=C(C)C1C(C)(C)C");
        this.fragmenter.setMolecule(mol1);
        Assert.assertTrue(compareLists(new List[]{List.of(2,3),List.of(6,7),List.of(11,12),List.of(14,15),List.of(14,16),List.of(6,8,9),List.of(13,14,17),List.of(4,2,1,0), new ArrayList()},
                this.fragmenter.getBranches()));
    }
    //</editor-fold>
    //<editor-fold desc="cutChains">
    //<editor-fold desc="Main Branch">

    /**
     *
     * @throws CDKException
     */
    @Test
    public void testCutChainsMainBranchWithoutRemainder () throws CDKException {
        IAtomContainer mol1 = this.sp.parseSmiles("CCCCCCCCCCCC");
        this.fragmenter.setMolecule(mol1);
        this.fragmenter.setFragmentationSettings(3,3,false);
        Assert.assertTrue(compareLists(new String[]{"CCC", "CCC", "CCC", "CCC"},
                getSmiles(this.fragmenter.getIAtomContainer())));
    }
    @Test
    public void testCutChainsMainBranchWithRemainderMinCut () throws CDKException {
        IAtomContainer mol1 = this.sp.parseSmiles("CCCCCCCCCCC");
        this.fragmenter.setMolecule(mol1);
        this.fragmenter.setFragmentationSettings(3,0,false);
        Assert.assertTrue(compareLists(new String[]{"CCC", "CCC", "C(CC)CC"},
                getSmiles(this.fragmenter.getIAtomContainer())));
    }
    @Test
    public void testCutChainsMainBranchWithRemainderMaxCut () throws CDKException {
        IAtomContainer mol1 = this.sp.parseSmiles("CCCCCCCCCCC");
        this.fragmenter.setMolecule(mol1);
        this.fragmenter.setFragmentationSettings(0,3,false);
        Assert.assertTrue(compareLists(new String[]{"CCC", "CCC", "CCC", "CC"},
                getSmiles(this.fragmenter.getIAtomContainer())));
    }
    @Test
    public void testCutChainsMainBranchMinCutMaxCut () throws CDKException {
        IAtomContainer mol1 = this.sp.parseSmiles("CCCCCCCCCCC");
        this.fragmenter.setMolecule(mol1);
        this.fragmenter.setFragmentationSettings(2,5,false);
        Assert.assertTrue(compareLists(new String[]{"CCCCC", "C(C)CCCC"},
                getSmiles(this.fragmenter.getIAtomContainer())));
    }
    @Test
    public void testCutChainsMainBranchMultipleBonds () throws CDKException {
        IAtomContainer mol1 = this.sp.parseSmiles("C=CC=CCC#CCC");
        this.fragmenter.setMolecule(mol1);
        this.fragmenter.setFragmentationSettings(2,3,false);
        Assert.assertTrue(compareLists(new String[]{"C#CCC", "C=CC", "C=C"},
                getSmiles(this.fragmenter.getIAtomContainer())));
    }
    @Test
    public void testCutChainsMainBranchAllenes () throws CDKException {
        IAtomContainer mol1 = this.sp.parseSmiles("CCC=C=CCCCC");
        this.fragmenter.setMolecule(mol1);
        this.fragmenter.setFragmentationSettings(2,3,false);
        Assert.assertTrue(compareLists(new String[]{"CCC", "C=C=CC", "CC"},
                getSmiles(this.fragmenter.getIAtomContainer())));
    }
    //</editor-fold>
    //<editor-fold desc="Sub Branch">
    /**
     *
     */
    @Test
    public void testCutChainsSubBranchDoubleBondAtBranching () throws CDKException {
        IAtomContainer mol1 = this.sp.parseSmiles("CCCCCCC(=CCCCC)CCCCCCC");
        this.fragmenter.setMolecule(mol1);
        this.fragmenter.setFragmentationSettings(2,3,false);
        Assert.assertTrue(compareLists(new String[]{"CCC", "CCC", "CCC", "CC", "C(C)CC", "CC(=C)C"},
                getSmiles(this.fragmenter.getIAtomContainer())));
    }
    @Test
    public void testCutChainsSubBranchAlleneAtBranching () throws CDKException {
        IAtomContainer mol1 = this.sp.parseSmiles("CCCCCCC(=C=CCCC)CCCCCCC");
        this.fragmenter.setMolecule(mol1);
        this.fragmenter.setFragmentationSettings(2,3,false);
        Assert.assertTrue(compareLists(new String[]{"CCC", "CCC", "CCC", "CCC", "CC", "CC(=C=C)C"},
                getSmiles(this.fragmenter.getIAtomContainer())));
    }
    @Test
    public void testCutChainsSubBranchIsPreservingTertiaryQuaternaryCarbons () throws CDKException {
        IAtomContainer mol1 = this.sp.parseSmiles("CCC(CC)(CCC)CC(C)CC");
        this.fragmenter.setMolecule(mol1);
        this.fragmenter.setFragmentationSettings(2,3,true);
        Assert.assertTrue(compareLists(new String[]{"CC(CC)(CC)CC(C)CC", "CC"},
                getSmiles(this.fragmenter.getIAtomContainer())));
    }
    @Test
    public void testCutChainsSubBranchNoMinCutMaxCut () throws CDKException {
        IAtomContainer mol1 = this.sp.parseSmiles("CCC(CC)(CCC)CC(C)CC");
        this.fragmenter.setMolecule(mol1);
        this.fragmenter.setFragmentationSettings(0,0,false);
        Assert.assertTrue(compareLists(new String[]{"CCCCCCCC", "CC", "CC", "C"},
                getSmiles(this.fragmenter.getIAtomContainer())));
    }
    @Test
    public void testCutChainsSubBranchNoMinCutMaxCutIsPreservingTertiaryQuaternaryCarbons () throws CDKException {
        IAtomContainer mol1 = this.sp.parseSmiles("CCC(CC)(CCC)CC(C)CC");
        this.fragmenter.setMolecule(mol1);
        this.fragmenter.setFragmentationSettings(0,0,true);
        Assert.assertTrue(compareLists(new String[]{"CCCC(C)(C)CC(C)CC", "C", "C"},
                getSmiles(this.fragmenter.getIAtomContainer())));
    }
    //</editor-fold>
    //</editor-fold>

    //<editor-fold desc="cutRings">
    /**
     *
     */
    @Test
    public void testCutRingsIsNotPreservingTertiaryQuaternaryCarbons () throws CDKException {
        IAtomContainer mol1 = this.sp.parseSmiles("C2CCC(CCCCCC1CCCCC1)CC2");
        this.fragmenter.setMolecule(mol1);
        this.fragmenter.setFragmentationSettings(2,3,false);
        Assert.assertTrue(compareLists(new String[]{"C1CCCCC1", "CCC", "CC", "C1CCCCC1"},
                getSmiles(this.fragmenter.getIAtomContainer())));
    }
    @Test
    public void testCutRingsIsPreservingTertiaryQuaternaryCarbons () throws CDKException {
        IAtomContainer mol1 = this.sp.parseSmiles("C2CCC(CCCCCC1CCCCC1)CC2");
        this.fragmenter.setMolecule(mol1);
        this.fragmenter.setFragmentationSettings(2,3,true);
        Assert.assertTrue(compareLists(new String[]{"C1CCC(C)CC1", "CC1CCCCC1", "CCC"},
                getSmiles(this.fragmenter.getIAtomContainer())));
    }
    //</editor-fold>
    //</editor-fold>
    //<editor-fold desc="Test Example Molecules">

    //</editor-fold>
}