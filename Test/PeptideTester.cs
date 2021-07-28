using NUnit.Framework;
using System;
using System.Collections.Generic;
using PTMStoichiometry;
using System.Linq;

namespace Test
{
    [TestFixture]
    class PeptideTester
    {
        
        private static readonly object[] _testFlashLFQGetter =
       {
            new object[] {

                @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-1000PeptidesProteinAlphabetized.txt",
                @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt", 
                3, 0,
                "AGLGTGAAGGIGAGRTRAPSLASSSGSDK", 
                "AGLGTGAAGGIGAGRTRAPSLASSSGSDK", 
                new List<string> { "A0A087WPF7" }, 
                "Auts2", 
                "Mus musculus", 
                18, 12,
                "Intensity_1DLC122419QE_ZD_Kidney_global_#4-calib", 
                "Global Control",
                3195479.843,
                true, true 
            },
            new object[] {

                @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-1000PeptidesProteinAlphabetized.txt",
                @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt",
                3, 1,
                "PSLDQIYTQFK",
                "PSLDQIYTQFK",
                new List<string> { "A0A0U1RPR8" },
                "Gucy2d",
                "Mus musculus",
                18, 6,
                "Intensity_1DLC122419QE_ZD_Kidney_global_#15-calib",
                "Global Day 7",
                2189433.856,
                true, false
            },
            new object[] {

                @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-1000PeptidesProteinAlphabetized.txt",
                @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt",
                3, 2,
                "ASVPSSEAGAWEVAASDIEPESRDRR",
                "ASVPSSEAGAWEVAASDIEPESRDRR",
                new List<string> { "A0A140LHF2" },
                "Vsig10l2",
                "Mus musculus",
                18, 2,
                "Intensity_1DLC122419QE_ZD_Kidney_global_#11-calib",
                "Global Day 2",
                0,
                true, false
            },
            new object[] {

                @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-1000PeptidesProteinAlphabetized.txt",
                @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt",
                3, 10,
                "LKHEC[Common Fixed:Carbamidomethyl on C]GAAFTSK",
                "LKHECGAAFTSK",
                new List<string> { "A2A432" },
                "Cul4b",
                "Mus musculus",
                18, 7,
                "Intensity_1DLC122419QE_ZD_Kidney_global_#16-calib",
                "Global Day 7",
                2088182.486,
                true, true
            }
        };

        /// <summary>
        /// Test checking that Peptide is reading in information from FlashLFQ lines and setting isUnique correctly
        /// </summary>
        /// <param name="filePath">file path to peptide data</param>
        /// <param name="groupPath">file path to group data</param>
        /// <param name="reqNumPepMeasurements">min num of peptide intensities that must be observed</param>
        /// <param name="lineToCompare">line to compare from peptide data</param>
        /// <param name="sequence">correct full sequence</param>
        /// <param name="baseSequence">correct base sequence</param>
        /// <param name="proteinGroup">correct protein group(s)</param>
        /// <param name="geneName">correct gene name</param>
        /// <param name="organism">correct organism</param>
        /// <param name="numberIntensities">correct number of intensities in the line</param>
        /// <param name="intensityColToCompare">column of the intensity to compare</param>
        /// <param name="fileName">correct file name of intensity</param>
        /// <param name="groupID">correct group ID of intensity</param>
        /// <param name="intensityVal">correct intensity value</param>
        /// <param name="isUnique">correct isUnique bool</param>
        /// <param name="detectedMinNum">correct detectedMinNum bool</param>
        [Test]
        [TestCaseSource("_testFlashLFQGetter")]
       
        public void Peptide_FlashLFQGetter_Pass(string filePath, string groupPath, int reqNumPepMeasurements, int lineToCompare, string sequence, string baseSequence,
            List<string> proteinGroup, string geneName, string organism, int numberIntensities, int intensityColToCompare, string fileName,
            string groupID, double intensityVal, bool isUnique, bool detectedMinNum)
        {
            //check that valid entries in constuctor lead to object in getter
            List<Peptide> pepsInFile = PeptideReader.ReadTsv(filePath, groupPath, reqNumPepMeasurements, 5, "FlashLFQ");
            Assert.AreEqual(sequence, pepsInFile[lineToCompare].Sequence);
            Assert.AreEqual(baseSequence, pepsInFile[lineToCompare].BaseSeq);
            Assert.AreEqual(proteinGroup, pepsInFile[lineToCompare].ProteinGroup);
            //Assert.AreEqual(proteinGroup, pepsInFile[lineToCompare].ProteinGroup);
            Assert.AreEqual(geneName, pepsInFile[lineToCompare].GeneName);
            Assert.AreEqual(numberIntensities, pepsInFile[lineToCompare].Intensities.Count);
            Assert.AreEqual(detectedMinNum, pepsInFile[lineToCompare].DetectedMinNum);

            //check setting isUnique
            //Extensions.IncludeSharedPeptides(pepsInFile, true);
            Assert.AreEqual(isUnique, pepsInFile[lineToCompare].IsUnique);

            //check intensities
            Assert.AreEqual(fileName, pepsInFile[lineToCompare].Intensities[intensityColToCompare].FileName);
            Assert.AreEqual(groupID, pepsInFile[lineToCompare].Intensities[intensityColToCompare].GroupID);
            Assert.That(intensityVal, Is.EqualTo(pepsInFile[lineToCompare].Intensities[intensityColToCompare].IntensityVal).Within(0.001));
        }

        /// <summary>
        /// Test which checks that intensitities are being read in correctly by reading in the same file with two different
        /// group files and checking that the intensities remain the same except for the group data  
        /// </summary>
        /// <param name="filePath">file path to peptide data</param>
        /// <param name="groupPath1">file path to first grouping data</param>
        /// <param name="groupPath2">file path to second grouping data</param>
        /// <param name="reqNumPepMeasurements">min num of peptide intensities that must be observed</param>
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-1000PeptidesProteinAlphabetized.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126TwoGroups.txt", 3)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-5000PeptidesProteinAlphabetized.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126TwoGroups.txt", 3)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-10000PeptidesProteinAlphabetized.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126TwoGroups.txt", 3)]
        public void Peptide_FlashLFQIntensity_Pass(string filePath, string groupPath1, string groupPath2, int reqNumPepMeasurements)
        {
            List<Peptide> pepsGroup1 = PeptideReader.ReadTsv(filePath, groupPath1, reqNumPepMeasurements, 5, "FlashLFQ");
            List<Peptide> pepsGroup2 = PeptideReader.ReadTsv(filePath, groupPath2, reqNumPepMeasurements, 5, "FlashLFQ");

            Assert.AreEqual(pepsGroup1.Count(), pepsGroup2.Count());
            Assert.AreEqual(pepsGroup1.Select(p => p.Sequence), pepsGroup2.Select(p => p.Sequence));

            for (int i = 0; i < pepsGroup1.Count(); i++)
            {
                Assert.AreEqual(pepsGroup1[i].Intensities.Count(), pepsGroup2[i].Intensities.Count());
                Assert.AreEqual(pepsGroup1[i].Intensities.Select(p => p.FileName), pepsGroup2[i].Intensities.Select(p => p.FileName));
                Assert.AreEqual(pepsGroup1[i].Intensities.Select(p => p.IntensityVal), pepsGroup2[i].Intensities.Select(p => p.IntensityVal));
                Assert.AreEqual(pepsGroup1[i].Intensities.Select(p => p.FileName), pepsGroup2[i].Intensities.Select(p => p.FileName));

            }
        }

        /// <summary>
        /// Test to check that DetectedMinNum is working correctly by adjusting reqNumPepMeasurements so different numbers of Peptides will qualify
        /// </summary>
        /// <param name="filePath">file path to the Peptides</param>
        /// <param name="groupPath">file path to the group data</param>
        /// <param name="reqNumPepMeasurements">min num of peptide intensities that must be observed</param>
        /// <param name="numPeptidesMeetMeasurementReq">correct number of Peptides that meet the requirement with the set DetectedMinNum</param>
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-1000PeptidesProteinAlphabetized.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt", 
            0, 1000)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-1000PeptidesProteinAlphabetized.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt",
            7, 0)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-1000PeptidesProteinAlphabetized.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt",
            6, 265)]
        public void Peptide_DetectedMinNum_Pass(string filePath, string groupPath, int reqNumPepMeasurements, int numPeptidesMeetMeasurementReq)
        {
            List<Peptide> pepsInFile = PeptideReader.ReadTsv(filePath, groupPath, reqNumPepMeasurements, 5, "FlashLFQ");
            Assert.AreEqual(numPeptidesMeetMeasurementReq, pepsInFile.Where(p => p.DetectedMinNum).Count());
        }
       
    }
}
