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
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.ArrayList;
import java.util.List;

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
    @Test
    public void testMethod() throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        testMinLength();
        testMaxLength();
    }

    public void testMinLength () throws InvalidSmilesException {
        this.sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer mol1 = sp.parseSmiles("CCCCCCCCCC");
        IAtomContainer mol2 = sp.parseSmiles("CCCC(CCCCC)CCCCC");
        this.fragmenter = new AlkylFragmenter();
        fragmenter.setMolecule(mol1);
        List<List<Integer>> actualIndices = this.fragmenter.chainCutter(mol1, this.fragmenter.branchCutter(mol1), 5, 0, false);
        List<List<Integer>> expectedIndices = List.of(List.of(5, 6, 7, 8, 9), List.of(0, 1, 2, 3, 4));
        Assert.assertEquals(expectedIndices, actualIndices);
        System.out.println(5);
        actualIndices = this.fragmenter.chainCutter(mol1, this.fragmenter.branchCutter(mol1), 4, -1, false);
        expectedIndices = List.of(List.of(0, 1, 2, 3), List.of(4, 5, 6, 7, 8, 9));
        Assert.assertEquals(expectedIndices, actualIndices);
        System.out.println(4);
        actualIndices = this.fragmenter.chainCutter(mol1, this.fragmenter.branchCutter(mol1), 3, -1, false);
        expectedIndices = List.of(List.of(3, 4, 5), List.of(0, 1, 2), List.of(6, 7, 8, 9));
        Assert.assertEquals(expectedIndices, actualIndices);
        System.out.println(3);
        System.out.println("---");
        this.fragmenter.setMolecule(mol2);
        actualIndices = this.fragmenter.chainCutter(mol2, this.fragmenter.branchCutter(mol2), 5, 0, false);
        expectedIndices = List.of(List.of(8, 7, 6, 5, 4), List.of(3, 9, 10, 11, 12, 2, 1, 0, 13));
        Assert.assertEquals(expectedIndices, actualIndices);
        System.out.println(5);
        actualIndices = this.fragmenter.chainCutter(mol2, this.fragmenter.branchCutter(mol2), 4, 0, false);
        expectedIndices = List.of(List.of(8, 7, 6, 5), List.of(4, 3, 9, 10, 2, 1, 0, 11, 12, 13));
        Assert.assertEquals(expectedIndices, actualIndices);
        System.out.println(4);
        actualIndices = this.fragmenter.chainCutter(mol2, this.fragmenter.branchCutter(mol2), 3, 0, false);
        expectedIndices = List.of(List.of(5, 4, 3), List.of(8, 7, 6), List.of(2, 1, 0), List.of(9, 10, 11, 12, 13));
        Assert.assertEquals(expectedIndices, actualIndices);
        System.out.println(3);
    }

    public void testMaxLength () throws InvalidSmilesException {
        this.sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer mol1 = sp.parseSmiles("CCCCCCCCCC");
        IAtomContainer mol2 = sp.parseSmiles("CCCC(CC)CCCCC");
        this.fragmenter = new AlkylFragmenter();
        fragmenter.setMolecule(mol1);
        List<List<Integer>> actualIndices = this.fragmenter.chainCutter(mol1, this.fragmenter.branchCutter(mol1), 1, 5, false);
        List<List<Integer>> expectedIndices = List.of(List.of(5,6,7,8,9), List.of(0,1,2,3,4));
        Assert.assertEquals(expectedIndices, actualIndices);
        System.out.println(5);
        actualIndices = this.fragmenter.chainCutter(mol1, this.fragmenter.branchCutter(mol1), 2, 4, false);
        expectedIndices = List.of(List.of(8, 9), List.of(4, 5, 6, 7), List.of(0, 1, 2, 3));
        Assert.assertEquals(expectedIndices, actualIndices);
        System.out.println(4);
        actualIndices = this.fragmenter.chainCutter(mol1, this.fragmenter.branchCutter(mol1), 2, 3, false);
        expectedIndices = List.of(List.of(9), List.of(6, 7, 8), List.of(3, 4, 5), List.of(0, 1, 2));
        Assert.assertEquals(expectedIndices, actualIndices);
        System.out.println(3);
        System.out.println("---");
        this.fragmenter.setMolecule(mol2);
        actualIndices = this.fragmenter.chainCutter(mol2, this.fragmenter.branchCutter(mol2), 1, 5, false);
        expectedIndices = List.of(List.of(7, 8, 9, 10), List.of(0, 1, 2, 3, 6), List.of(4, 5));
        Assert.assertEquals(expectedIndices, actualIndices);
        System.out.println(5);
        actualIndices = this.fragmenter.chainCutter(mol2, this.fragmenter.branchCutter(mol2), 1, 4, false);
        expectedIndices = List.of(List.of(10), List.of(6, 7, 8, 9), List.of(0, 1, 2, 3), List.of(4, 5));
        Assert.assertEquals(expectedIndices, actualIndices);
        System.out.println(4);
        actualIndices = this.fragmenter.chainCutter(mol2, this.fragmenter.branchCutter(mol2), 1, 3, false);
        expectedIndices = List.of(List.of(8, 9, 10), List.of(3, 6, 7), List.of(0, 1, 2), List.of(4, 5));
        Assert.assertEquals(expectedIndices, actualIndices);
        System.out.println(3);
    }
}