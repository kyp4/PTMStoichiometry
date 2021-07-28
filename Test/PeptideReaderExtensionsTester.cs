using NUnit.Framework;
using System;
using System.Collections.Generic;
using PTMStoichiometry;
using System.Linq;
using System.IO;
using System.Text;

namespace Test
{
    [TestFixture]
    class PeptideReaderExtensionsTester
    {

        /// <summary>
        /// Test that checks that ReadTsv returns the correct number of Peptides
        /// </summary>
        /// <param name="filePath">file path to the file with the peptide data</param>
        /// <param name="groupPath">file path to the file with the group data</param>
        /// <param name="reqNumPepMeasurements">min num of peptide intensities that must be observed (non zero -> MS or MSMS detection)</param>
        /// <param name="peptideCount">the number of Peptides in the file</param>
        [Test]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-1000PeptidesProteinAlphabetized.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt", 3, 1000)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-5000PeptidesProteinAlphabetized.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt", 3, 5000)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-10000PeptidesProteinAlphabetized.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt", 3, 10000)]
        public void PeptideReader_ReadTsv_Pass(string filePath, string groupPath, int reqNumPepMeasurements, int peptideCount)
        {
            //check throws error if file references no good
            //var noFilePeptideReaderTest = PeptideReader.ReadTsv("", filePathGroupsOnePeptideTest);
            //Assert.Throws(, onePeptideReaderTest.Count);

            List<Peptide> pepsInFile = PeptideReader.ReadTsv(filePath, groupPath, reqNumPepMeasurements, 5, "FlashLFQ");
            Assert.AreEqual(peptideCount, pepsInFile.Count());

            //test group dictonaries & lists!
        }

        [Test]
        [TestCase(1, 1, 2, true, 3, 3, 0.75, false, 3, 3, true, 0.05, null, "FlashLFQ",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides.tsv",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\pvals\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-PeptideAnalysis-Correctedpvals.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\pvals\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-PTMAnalysis-Correctedpvals.txt")]
        public void BenjaminiHochbergTester(
            int reqNumUnmodPeptides,
            int reqNumModPeptides,
            int reqNumOfPepeptides,
            bool useBaselinePeptides,
            int reqNumBaselinePeptides,
            int reqNumBaselineMeasurements,
            double correlationCutOff,
            bool compareUnmod,
            int minNumStoichiometries,
            int reqNumPepMeasurements,
            bool groupPepsForPValCalc,
            double alpha,
            string groupToCompare,
            string dataType,

            string filepathpeptides,
            string filepathgroups,

            string filepathcorrectpvalsprot,
            string filepathcorrectpvalsptm)
        {

            if (File.ReadAllLines(filepathpeptides, Encoding.UTF8)[0].Split("\t")[4] == "Organism")
            {
                dataType = "FlashLFQ";
            }
            else
            {
                dataType = "MaxQuant";
            }
            string subdirectory = "PXD003881-20210722a";
            string peptidestoichiometryfileout = subdirectory + "PeptideAnalysis";
            string ptmstoichiometryfileout = subdirectory + "PTMAnalysis";
            string paramsfile = subdirectory + "params";

            int intensityIndex = 5;
            if (dataType == "MaxQuant")
            {
                intensityIndex = PeptideReader.IndexFind(filepathpeptides, "Intensity");
            }


            Dictionary<string, string> groups = PeptideReader.GetGroups(filepathgroups);
            List<string> groupsList = PeptideReader.GetGroupList(filepathgroups);
            List<Peptide> testPeptide = PeptideReader.ReadTsv(filepathpeptides, filepathgroups, reqNumPepMeasurements, intensityIndex, dataType);
            
            var proteins = testPeptide.Where(p => p.ProteinGroup.Count() == 1).Select(p => p.ProteinGroup).Distinct().ToList();

            List<string> proteinList = new List<string>();
            foreach (List<string> prot in proteins)
            {
                proteinList.Add(prot[0]);
            }
            proteinList = proteinList.Distinct().ToList();
            List<ProteinGroup> testProt = new List<ProteinGroup>();

            if (groupToCompare != null)
            {
                for (int i = 0; i < proteinList.Count(); i++)
                {
                    testProt.Add(new ProteinGroup(proteinList[i], testPeptide, groupsList, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides,
                    useBaselinePeptides, reqNumBaselinePeptides, reqNumBaselineMeasurements, correlationCutOff, compareUnmod, minNumStoichiometries, groupToCompare));
                }
            }
            else
            {
                for (int i = 0; i < proteinList.Count(); i++)
                {
                    testProt.Add(new ProteinGroup(proteinList[i], testPeptide, groupsList, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides,
                    useBaselinePeptides, reqNumBaselinePeptides, reqNumBaselineMeasurements, correlationCutOff, compareUnmod, minNumStoichiometries));
                }
            }

            testProt = testProt.Where(p => p.ProteinPairwiseComparisons != null).Distinct().ToList(); //hmmm
            Extensions.CalcCorrectedPValue(testProt, groupPepsForPValCalc, alpha);
            List<ProteinGroup> ProteinsToUse = testProt.Where(p => p.useProt).ToList();

            List<double> pvalsPTM = new List<double>();
            List<double> pvalsProt = new List<double>();
            foreach (ProteinGroup prot in ProteinsToUse)
            {
               foreach (PairwiseCompairison pc in prot.ProteinPairwiseComparisons)
               {
                    pvalsProt.Add(pc.CorrectedpVal);
               }
               foreach (PairwiseCompairison ptmc in prot.PTMPairwiseCompairisons)
               {
                    pvalsPTM.Add(ptmc.CorrectedpVal);
               }
            }

            pvalsProt.Sort();
            pvalsPTM.Sort();

            List<string> correctpvalsProtS = File.ReadLines(filepathcorrectpvalsprot).ToList();
            List<double> correctpvalsProt = correctpvalsProtS.Select(p => double.Parse(p)).ToList();
            correctpvalsProt.Sort();

            List<string> correctpvalsPTMS = File.ReadLines(filepathcorrectpvalsptm).ToList();
            List<double> correctpvalsPTM = correctpvalsPTMS.Select(p => double.Parse(p)).ToList();
            correctpvalsPTM.Sort();

            Assert.That(pvalsProt, Is.EqualTo(correctpvalsProt).Within(0.1));
            Assert.That(pvalsPTM, Is.EqualTo(correctpvalsPTM).Within(0.1));
        }

    }
}
