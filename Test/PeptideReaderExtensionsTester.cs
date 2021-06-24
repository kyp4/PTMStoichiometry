using NUnit.Framework;
using System;
using System.Collections.Generic;
using PTMStoichiometry;
using System.Linq;
using System.IO;

namespace Test
{
    [TestFixture]
    class PeptideReaderExtensionsTester
    {

        //test that PeptideReader reads in the correct number of Peptides
        [Test]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPosphoSmallTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 31)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPosphoLargeTest.tsv",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 39869)]
        public void PeptideReader_ReadTsv_Pass(string filePath, string groupPath, int reqNumPepMeasurements, int peptideCount)
        {
            //check throws error if file references no good
            //var noFilePeptideReaderTest = PeptideReader.ReadTsv("", filePathGroupsOnePeptideTest);
            //Assert.Throws(, onePeptideReaderTest.Count);

            List<Peptide> pepsInFile = PeptideReader.ReadTsv(filePath, groupPath, reqNumPepMeasurements, 5, "FlashLFQ");
            Assert.AreEqual(peptideCount, pepsInFile.Count());

            //test group dictonaries & lists!
        }

        //test that all shared peptides are being identified by IncludeSharedPeptidesTester
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinySharedPeptideTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3,
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt", 2)]
        [Test]
        public void Extensions_IncludeSharedPeptides_Pass(string filePathSharedPeptides, string groupPath, int reqNumPepMeasurements, string filePathUniquePeptides, int numRazorPep)
        {
            List<Peptide> sharedPeps = PeptideReader.ReadTsv(filePathSharedPeptides, groupPath, reqNumPepMeasurements, 5, "FlashLFQ");
            List<Peptide> uniquePeps = PeptideReader.ReadTsv(filePathUniquePeptides, groupPath, reqNumPepMeasurements, 5, "FlashLFQ");

            //include razor peptides
            //sharedPeps = Extensions.IncludeSharedPeptides(sharedPeps, true);
            Assert.AreNotEqual(uniquePeps.Count(), sharedPeps.Count());
            //remove razor peptides
            //sharedPeps = Extensions.IncludeSharedPeptides(sharedPeps, false);

            Assert.AreEqual(uniquePeps.Count(), sharedPeps.Count());
        }

        [Test]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPosphoSmallTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 1, 4, 5, false, 0.05,
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\PTMStoichiometrySmallPhosphoData-pvals-ungroupedR.txt", 
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\PTMStoichiometrySmallPhosphoData-pvals-groupedR.txt")]
        public void BenjaminiHochbergTester(string filepathpeptides, string filepathgroups,
            int reqNumUnmodPeptides, int reqNumModPeptides, int reqNumOfPepeptides, int reqNumBaselineMeasurements, 
            Boolean useBaselinePeptides, double alpha,
            string correctPValsUngroupedFile, string correctPValsGroupedFile)
        {
            List<string> correctPValsUngroupedS = File.ReadLines(correctPValsUngroupedFile).ToList();
            List<double> correctPValsUngrouped = correctPValsUngroupedS.Select(p => double.Parse(p)).ToList();
            List<string> correctPValsGroupedS = File.ReadLines(correctPValsGroupedFile).ToList();
            List<double> correctPValsGrouped = correctPValsGroupedS.Select(p => double.Parse(p)).ToList();
            
            
            Dictionary<string, string> groups = PeptideReader.GetGroups(filepathgroups);
            List<string> groupsList = PeptideReader.GetGroupList(filepathgroups);
            List<Peptide> testPeptide = PeptideReader.ReadTsv(filepathpeptides, filepathgroups, 3, 5, "FlashLFQ");
            //testPeptide = Extensions.IncludeSharedPeptides(testPeptide, false);

            //group peptides by protein
            var proteins = testPeptide.Select(p => p.ProteinGroup).Distinct().ToArray();
            List<ProteinGroup> testProt = new List<ProteinGroup>();
            for (int i = 0; i < proteins.Length; i++)
            {
                testProt.Add(new ProteinGroup(proteins[i][0], testPeptide, groupsList, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides,
                    useBaselinePeptides, reqNumUnmodPeptides, reqNumBaselineMeasurements, 0.5, false, 3));
            }
            var correctPValsUnGroupedProt = testProt.Where(p => p.ProteinPairwiseComparisons != null).ToList();
            var correctPValsGroupedProt = testProt.Where(p => p.ProteinPairwiseComparisons != null).ToList();
            Extensions.CalcCorrectedPValue(correctPValsUnGroupedProt, false, alpha);
            Extensions.CalcCorrectedPValue(correctPValsGroupedProt, true, alpha);

            List<double> pValsGroupedPeps = new List<double>();
            List<double> pValsUngroupedPeps = new List<double>();
            for (int i=0; i < correctPValsUnGroupedProt.Count(); ++i)
            {
                pValsGroupedPeps.AddRange(correctPValsGroupedProt[i].ProteinPairwiseComparisons.Select(p => p.CorrectedpVal));
                pValsUngroupedPeps.AddRange(correctPValsUnGroupedProt[i].ProteinPairwiseComparisons.Select(p => p.CorrectedpVal));
            }
            pValsGroupedPeps.Sort();
            correctPValsGrouped.Sort();
            Assert.That(pValsUngroupedPeps, Is.EqualTo(correctPValsUngrouped).Within(0.001));
            Assert.That(pValsGroupedPeps, Is.EqualTo(correctPValsGrouped).Within(0.001));
        }

    }
}
