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

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.xmlcml.euclid.Int;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class ConsoleApplication {
    public static void main(String args[]) throws CDKException, IOException, CloneNotSupportedException {
        BufferedReader tmpBufferedReader = new BufferedReader(new InputStreamReader(System.in));
        AlkylFragmenter tmpFragmenter = new AlkylFragmenter();
        SmilesParser tmpSmilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        int tmpMinCut = 1;
        int tmpMaxCut = 3;
        String tmpMoleculeSmiles = "";
        String tmpPreferences = "";
        boolean tmpIsPreservingTertiaryQuaternaryCarbonAtoms = true;
        String tmpTask = "";
        while (!tmpTask.equals("4")) {
            System.out.println("What do you want to do?");
            System.out.println("1 Enter a Molecule SMILES String");
            System.out.println("2 Set Preferences");
            System.out.println("3 Import CSV File");
            System.out.println("4 Exit");
            tmpTask = tmpBufferedReader.readLine();
            if (tmpTask.equals("1") || tmpMoleculeSmiles.equals("") && !tmpTask.equals("3")) {
                System.out.print("Molecule Smiles: ");
                tmpMoleculeSmiles = tmpBufferedReader.readLine();
                IAtomContainer tmpMol = tmpSmilesParser.parseSmiles(tmpMoleculeSmiles);
                tmpFragmenter.setMolecule(tmpMol);
            }
            while (tmpTask.equals("2") && !tmpPreferences.equals("4")) {
                System.out.println("Which preference do you want to change?");
                System.out.println("1 Minimum Chain Length (current: " + Integer.toString(tmpMinCut)+")");
                System.out.println("2 Maximum Chain Length (current: " + Integer.toString(tmpMaxCut)+")");
                System.out.println("3 Preserve Tertiary and Quaternary Carbon Atoms (current: " + Boolean.toString(tmpIsPreservingTertiaryQuaternaryCarbonAtoms)+")");
                System.out.println("4 Back");
                tmpPreferences = tmpBufferedReader.readLine();
                if (tmpPreferences.equals("1")) {
                    System.out.print("New Minimum Chain Length: ");
                    tmpMinCut = Integer.parseInt(tmpBufferedReader.readLine());
                }
                if (tmpPreferences.equals("2")) {
                    System.out.print("New Maximum Chain Length: ");
                    tmpMaxCut = Integer.parseInt(tmpBufferedReader.readLine());
                }
                if (tmpPreferences.equals("3")) {
                    System.out.println("Preserve Tertiary and Quaternary Carbon Atoms?\n1 Yes\n2 No");
                    String tmpPreserveTertiaryQuaternaryCarbonOption = tmpBufferedReader.readLine();
                    if (tmpPreserveTertiaryQuaternaryCarbonOption.equals("1")) {
                        tmpIsPreservingTertiaryQuaternaryCarbonAtoms = true;
                    }
                    if (tmpPreserveTertiaryQuaternaryCarbonOption.equals("2")) {
                        tmpIsPreservingTertiaryQuaternaryCarbonAtoms = false;
                    }
                }
            }
            tmpPreferences = "";
            tmpFragmenter.fragmentationSettings(tmpMinCut, tmpMaxCut, tmpIsPreservingTertiaryQuaternaryCarbonAtoms);
            System.out.println(tmpFragmenter.getSmiles());
        }
    }
}