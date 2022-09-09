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
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 */
public class AlkylFragmenterTest {
    /**
     *
     * @throws Exception
     */

    SmilesParser sp;
    AlkylFragmenter fragmenter;
    Logger logger;
    public AlkylFragmenterTest () {
        sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        this.fragmenter = new AlkylFragmenter();
        logger = Logger.getLogger(AlkylFragmenterTest.class.getName());
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

    private boolean compareSmilesList (String[] aSmilesListExpected, List<String> aSmilesListActual) {
        List<String> tmpSmilesList1 = aSmilesListActual;
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

    private void testNumberOfAtoms () {
        int tmpNumberOfAtomsInFragmentList = 0;
        for (IAtomContainer tmpFragment : this.fragmenter.getIAtomContainer()) {
            tmpNumberOfAtomsInFragmentList += tmpFragment.getAtomCount();
        }
        Assert.assertEquals(this.fragmenter.getMolecule().getAtomCount(), tmpNumberOfAtomsInFragmentList);
    }
    @Test
    public void testPreserveCarbons () throws CDKException, CloneNotSupportedException {
        IAtomContainer mol1 = this.sp.parseSmiles("CCCC(C(C)(C)C)CCC(CC(C)(C(C)C))CC(CC)C(C)C");
        this.fragmenter.setMolecule(mol1);
        this.fragmenter.fragmentationSettings(1, 0, true);
        testNumberOfAtoms();
        Assert.assertTrue(compareSmilesList(new String[]{"C", "C", "C"}, this.fragmenter.getSmiles()));
    }

    @Test
    public void testMinLength () throws CDKException, CloneNotSupportedException {
        IAtomContainer mol1 = this.sp.parseSmiles("CCCCCCCCCC");
        IAtomContainer mol2 = this.sp.parseSmiles("CCCC(CCCCC)CCCCC");
        fragmenter.setMolecule(mol1);
        this.fragmenter.fragmentationSettings(5,0,false);
        testNumberOfAtoms();
        Assert.assertTrue(compareSmilesList(new String[]{"C", "C", "C"}, this.fragmenter.getSmiles()));
        this.fragmenter.fragmentationSettings(4,0,false);
        testNumberOfAtoms();
        Assert.assertTrue(compareSmilesList(new String[]{"C", "C", "C"}, this.fragmenter.getSmiles()));
        this.fragmenter.fragmentationSettings(3,0,false);
        testNumberOfAtoms();
        Assert.assertTrue(compareSmilesList(new String[]{"C", "C", "C"}, this.fragmenter.getSmiles()));
        this.fragmenter.setMolecule(mol2);
        this.fragmenter.fragmentationSettings(5,0,false);
        testNumberOfAtoms();
        Assert.assertTrue(compareSmilesList(new String[]{"C", "C", "C"}, this.fragmenter.getSmiles()));
        this.fragmenter.fragmentationSettings(4,0,false);
        testNumberOfAtoms();
        Assert.assertTrue(compareSmilesList(new String[]{"C", "C", "C"}, this.fragmenter.getSmiles()));
        this.fragmenter.fragmentationSettings(3,0,false);
        testNumberOfAtoms();
        Assert.assertTrue(compareSmilesList(new String[]{"C", "C", "C"}, this.fragmenter.getSmiles()));
    }
    @Test
    public void testMaxLength () throws CDKException, CloneNotSupportedException {
        IAtomContainer mol1 = this.sp.parseSmiles("CCCCCCCCCC");
        IAtomContainer mol2 = this.sp.parseSmiles("CCCC(CC)CCCCC");
        fragmenter.setMolecule(mol1);
        this.fragmenter.fragmentationSettings(1,5,false);
        testNumberOfAtoms();
        Assert.assertTrue(compareSmilesList(new String[]{"C", "C", "C"}, this.fragmenter.getSmiles()));
        this.fragmenter.fragmentationSettings(1,4,false);
        testNumberOfAtoms();
        Assert.assertTrue(compareSmilesList(new String[]{"C", "C", "C"}, this.fragmenter.getSmiles()));
        this.fragmenter.fragmentationSettings(1,3,false);
        testNumberOfAtoms();
        Assert.assertTrue(compareSmilesList(new String[]{"C", "C", "C"}, this.fragmenter.getSmiles()));
        this.fragmenter.setMolecule(mol2);
        this.fragmenter.fragmentationSettings(1,5,false);
        testNumberOfAtoms();
        Assert.assertTrue(compareSmilesList(new String[]{"C", "C", "C"}, this.fragmenter.getSmiles()));
        this.fragmenter.fragmentationSettings(1,4,false);
        testNumberOfAtoms();
        Assert.assertTrue(compareSmilesList(new String[]{"C", "C", "C"}, this.fragmenter.getSmiles()));
        this.fragmenter.fragmentationSettings(1,3,false);
        testNumberOfAtoms();
        Assert.assertTrue(compareSmilesList(new String[]{"C", "C", "C"}, this.fragmenter.getSmiles()));
    }
}