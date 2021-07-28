using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace PTMStoichiometry
{
    //class to write output
    public class WriteFile
    {
        /// <summary>
        /// Function to write peptide stoichiometry data to tsv 
        /// </summary>
        /// <param name="proteinGroups">list of protein groups</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        /// <param name="outputFolder">folder to write data to</param>
        /// <param name="fileName">name of file to write to</param>
        public static void StoichiometryPeptideDataWriter(List<ProteinGroup> proteinGroups, int minNumStoichiometries, string outputFolder, string fileName)
        {
            if (proteinGroups.Count == 0)
            { return; }
            var writtenFile = Path.Combine(outputFolder, fileName + ".txt");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                
                    output.WriteLine("Peptide\tModifications\tBaseline\tGroup1\tGroup2\tProtein\tTotalNumPeptidesProt\tMWStat\tMWPVal\tBHPVal\t" 
                        + "MedianStoichiometryGroup1\tMedianStoichiometryGroup2\tMinStoichiometryGroup1\tMinStoichiometryGroup2\t"
                        + "MaxStoichiometryGroup1\tMaxStoichiometryGroup2\tStoichiometriesGroup1\tStoichiometriesGroup2");
                    foreach (ProteinGroup prot in proteinGroups)
                    {
                        foreach (PairwiseCompairison pc in prot.ProteinPairwiseComparisons)
                        {
                            if (pc.PeptideStoichiometriesGroupOne.Count() > minNumStoichiometries && 
                                pc.PeptideStoichiometriesGroupTwo.Count() > minNumStoichiometries)
                            {
                                output.WriteLine(
                                    pc.PeptideOne.Sequence
                                    + "\t" + String.Join(";", pc.PeptideOne.PostTranslationalModifications.Select(p => p.Modification))
                                    + "\t" + String.Join(";", prot.BaselinePeptides.Select(p => p.Sequence))
                                    + "\t" + pc.GroupOne
                                    + "\t" + pc.GroupTwo
                                    + "\t" + String.Join(";", prot.ProteinName)
                                    + "\t" + prot.NumPeptidesInProtein.ToString()
                                    + "\t" + pc.MWStat.ToString()
                                    + "\t" + pc.MWPVal.ToString()
                                    + "\t" + pc.PeptideStoichiometriesGroupOneMedian.ToString()
                                    + "\t" + pc.PeptideStoichiometriesGroupTwoMedian.ToString()
                                    + "\t" + pc.PeptideStoichiometriesGroupOneMin.ToString()
                                    + "\t" + pc.PeptideStoichiometriesGroupTwoMin.ToString()
                                    + "\t" + pc.PeptideStoichiometriesGroupOneMax.ToString()
                                    + "\t" + pc.PeptideStoichiometriesGroupTwoMax.ToString()
                                    + "\t" + String.Join(";", pc.PeptideStoichiometriesGroupOne.Select(p => p.StoichiometryVal))
                                    + "\t" + String.Join(";", pc.PeptideStoichiometriesGroupTwo.Select(p => p.StoichiometryVal))
                                );
                            }
                        }
                    }
                
                
                
            }
        }

        /// <summary>
        /// Function to write ptm stoichiometry data to tsv 
        /// </summary>
        /// <param name="proteinGroups">list of protein groups</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        /// <param name="outputFolder">folder to write data to</param>
        /// <param name="fileName">name of file to write to</param>
        public static void StoichiometryPTMDataWriter(List<ProteinGroup> proteinGroups, int minNumStoichiometries, string outputFolder, string fileName)
        {
            if (proteinGroups.Count == 0)
            { return; }
            var writtenFile = Path.Combine(outputFolder, fileName + ".txt");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                
                    output.WriteLine("Modification\tPeptides\tProtein\tBaseline\tGroup1\tGroup2\tMWStat\tMWPVal\tBHPVal\t"
                        + "MedianStoichiometryGroup1\tMedianStoichiometryGroup2\tMinStoichiometryGroup1\tMinStoichiometryGroup2\t"
                        + "MaxStoichiometryGroup1\tMaxStoichiometryGroup2\tStoichiometriesGroup1\tStoichiometriesGroup2");
                    foreach (ProteinGroup prot in proteinGroups)
                    {
                        foreach (PairwiseCompairison pc in prot.PTMPairwiseCompairisons)
                        {
                            if (pc.PeptideStoichiometriesGroupOne.Count() > minNumStoichiometries &&
                                pc.PeptideStoichiometriesGroupTwo.Count() > minNumStoichiometries)
                            {
                                output.WriteLine(
                                    pc.PTM
                                    + "\t" + String.Join(";", pc.PeptidesWithPTM.Select(p => p.Sequence))
                                    + "\t" + String.Join(";", prot.ProteinName)
                                    + "\t" + String.Join(";", prot.BaselinePeptides.Select(p => p.Sequence))
                                    + "\t" + pc.GroupOne
                                    + "\t" + pc.GroupTwo
                                    + "\t" + pc.MWStat.ToString()
                                    + "\t" + pc.MWPVal.ToString()
                                    + "\t" + pc.PeptideStoichiometriesGroupOneMedian.ToString()
                                    + "\t" + pc.PeptideStoichiometriesGroupTwoMedian.ToString()
                                    + "\t" + pc.PeptideStoichiometriesGroupOneMin.ToString()
                                    + "\t" + pc.PeptideStoichiometriesGroupTwoMin.ToString()
                                    + "\t" + pc.PeptideStoichiometriesGroupOneMax.ToString()
                                    + "\t" + pc.PeptideStoichiometriesGroupTwoMax.ToString()
                                    + "\t" + String.Join(";", pc.PeptideStoichiometriesGroupOne.Select(p => p.StoichiometryVal))
                                    + "\t" + String.Join(";", pc.PeptideStoichiometriesGroupTwo.Select(p => p.StoichiometryVal))
                                );
                            }
                        }
                    }
            }
        }


        /// <summary>
        /// Function to write parameters for run to tsv 
        /// </summary>
        /// <param name="paramsfile">file name for paramter file</param>
        /// <param name="filepathpeptides">file path for peptide data</param>
        /// <param name="filepathgroups">file path for group data</param>
        /// <param name="directory">directory for output</param>
        /// <param name="stoichiometryfileout">where the stoichiometry outputs were written and file name</param>
        /// <param name="reqNumUnmodPeptides">min num of unmodified peptides that must be observed for a protein in order to consider the protein group (default=1)</param>
        /// <param name="reqNumModPeptides">min num of modified peptides that must be observed for a protein in order to consider the protein group (default=3)</param>
        /// <param name="reqNumOfPepeptides">min num of peptides that must be observed for a protein in order to consider the protein group </param>
        /// <param name="reqNumBaselinePeptides">min num of baseline peptides that must be observed for a protein in order to consider it (default=3)</param>
        /// <param name="reqNumBaselineMeasurements">min num of intensities (non zero -> MS or MSMS detection) that must be observed for in a 
        ///baseline peptide (non zero -> MS or MSMS detection) - increasing this value will decrease the number of baseline peptides  
        ///that are not observed in samples and therefore the number of non numeric stoichiometry values found in baseline case</param>
        /// <param name="correlationCutOff">min value at which two peptides will be considered to be correlated</param>
        /// <param name="compareUnmod">if false (default) only compare modified peptides to baseline, not unmodified peptides</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        /// <param name="reqNumPepMeasurements">min num of peptide intensities that must be observed (non zero -> MS or MSMS detection)</param>
        /// <param name="groupToCompare">single group to compare all other groups to (not a required parameter)</param>
        public static void ParamsWriter(string paramsfile, string filepathpeptides, string filepathgroups, string directory, string stoichiometryfileout,
            int reqNumUnmodPeptides = 1, int reqNumModPeptides = 3, int reqNumOfPepeptides = 4,
            int reqNumBaselinePeptides = 3, int reqNumBaselineMeasurements = 3, double correlationCutOff = 0.5,
            bool compareUnmod = false, int minNumStoichiometries = 3, int reqNumPepMeasurements = 3, string groupToCompare = null)
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
                    + "\n" + "reqNumBaselinePeptides:" + "\t" + reqNumBaselinePeptides
                    + "\n" + "reqNumBaselineMeasurements:" + "\t" + reqNumBaselineMeasurements
                    + "\n" + "correlationCutOff:" + "\t" + correlationCutOff
                    + "\n" + "compareUnmod:" + "\t" + compareUnmod
                    + "\n" + "minNumStoichiometries:" + "\t" + minNumStoichiometries
                    + "\n" + "reqNumPepMeasurements:" + "\t" + reqNumPepMeasurements
                    + "\n" + "groupToCompare:" + "\t" + groupToCompare
                    );

        }
    }
}
