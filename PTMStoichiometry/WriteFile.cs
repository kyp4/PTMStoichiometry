﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace PTMStoichiometry
{
    public class WriteFile
    {

        public static void StoichiometryDataWriter(List<ProteinGroup> proteinGroups, bool useBaselinePeptides, int minNumStoichiometries, string outputFolder, string fileName)
        {
            if (proteinGroups.Count == 0)
            { return; }
            var writtenFile = Path.Combine(outputFolder, fileName + ".txt");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                if (useBaselinePeptides)
                {
                    output.WriteLine("Peptide\tBaseline\tGroup1\tGroup2\tProtein\tTotalNumPeptidesProt\tMWStat\tMWPVal\tBHPVal\tStoichiometriesGroup1\tStoichiometriesGroup2");
                    foreach (ProteinGroup prot in proteinGroups)
                    {
                        foreach (PairwiseCompairison pc in prot.ProteinPairwiseComparisons)
                        {
                            if (pc.PeptideStoichiometriesGroupOne.Count() > minNumStoichiometries && 
                                pc.PeptideStoichiometriesGroupTwo.Count() > minNumStoichiometries)
                            {
                                output.WriteLine(
                                    pc.PeptideOne.Sequence
                                    + "\t" + String.Join(";", prot.BaselinePeptides.Select(p => p.Sequence))
                                    + "\t" + pc.GroupOne
                                    + "\t" + pc.GroupTwo
                                    + "\t" + prot.ProteinName
                                    + "\t" + prot.NumPeptidesInProtein.ToString()
                                    + "\t" + pc.MWStat.ToString()
                                    + "\t" + pc.MWPVal.ToString()
                                    + "\t" + pc.CorrectedpVal.ToString()
                                    + "\t" + String.Join(";", pc.PeptideStoichiometriesGroupOne.Select(p => p.StoichiometryVal))
                                    + "\t" + String.Join(";", pc.PeptideStoichiometriesGroupTwo.Select(p => p.StoichiometryVal))
                                );
                            }
                        }
                    }
                }
                else
                {
                    output.WriteLine("Peptide1\tPeptide2\tGroup1\tGroup2\tProtein\tTotalNumPeptidesProt\tMWStat\tMWPVal\tBHPVal\tStoichiometriesGroup1\tStoichiometriesGroup2");
                    foreach (ProteinGroup prot in proteinGroups)
                    {
                        foreach (PairwiseCompairison pc in prot.ProteinPairwiseComparisons)
                        {
                            if (pc.PeptideStoichiometriesGroupOne.Count() > minNumStoichiometries &&
                                pc.PeptideStoichiometriesGroupTwo.Count() > minNumStoichiometries)
                            {
                                output.WriteLine(
                                    pc.PeptideOne.Sequence
                                    + "\t" + pc.PeptideTwo.Sequence
                                    + "\t" + pc.GroupOne
                                    + "\t" + pc.GroupTwo
                                    + "\t" + prot.ProteinName
                                    + "\t" + prot.NumPeptidesInProtein.ToString()
                                    + "\t" + pc.MWStat.ToString()
                                    + "\t" + pc.MWPVal.ToString()
                                    + "\t" + pc.CorrectedpVal.ToString()
                                    + "\t" + String.Join(";", pc.PeptideStoichiometriesGroupOne.Select(p => p.StoichiometryVal))
                                    + "\t" + String.Join(";", pc.PeptideStoichiometriesGroupTwo.Select(p => p.StoichiometryVal))
                                );
                            }
                        }
                    }
                }
                
            }
        }

     

        public static void ParamsWriter(string paramsfile, string filepathpeptides, string filepathgroups, string directory, string stoichiometryfileout, 
            int reqNumUnmodPeptides = 1, int reqNumModPeptides = 3, int reqNumOfPepeptides = 4, bool useBaselinePeptides = true,
            int reqNumBaselinePeptides = 3, int reqNumBaselineMeasurements = 3, double correlationCutOff = 0.5, 
            bool compareUnmod = false, int minNumStoichiometries = 3, int reqNumPepMeasurements = 3, bool groupPepsForPValCalc = true, 
            double alpha = 0.05, string groupToCompare = null)
        {
            if (filepathpeptides == null || filepathgroups == null)
            { return; }
            var writtenFile = Path.Combine(directory, paramsfile + ".txt");
            using (StreamWriter output = new StreamWriter(writtenFile))

            output.WriteLine("Peptides Data:" + "\t" + filepathpeptides
                + "\n" + "Groups Data:" + "\t" + filepathgroups
                + "\n" + "Output Data:" + "\t" + directory + stoichiometryfileout
                + "\n" + DateTime.Now.ToString()
                + "\n" + "------------------------------------"
                + "\n" + "reqNumUnmodPeptides:" + "\t" + reqNumUnmodPeptides
                + "\n" + "reqNumModPeptides:" + "\t" + reqNumModPeptides
                + "\n" + "reqNumOfPepeptides:" + "\t" + reqNumOfPepeptides
                + "\n" + "useBaselinePeptides:" + "\t" + useBaselinePeptides
                + "\n" + "reqNumBaselinePeptides:" + "\t" + reqNumBaselinePeptides
                + "\n" + "reqNumBaselineMeasurements:" + "\t" + reqNumBaselineMeasurements
                + "\n" + "correlationCutOff:" + "\t" + correlationCutOff
                + "\n" + "compareUnmod:" + "\t" + compareUnmod
                + "\n" + "minNumStoichiometries:" + "\t" + minNumStoichiometries
                + "\n" + "reqNumPepMeasurements:" + "\t" + reqNumPepMeasurements
                + "\n" + "groupPepsForPValCalc:" + "\t" + groupPepsForPValCalc
                + "\n" + "alpha:" + "\t" + alpha
                + "\n" + "groupToCompare:" + "\t" + groupToCompare
                );
            



        }

    }
}