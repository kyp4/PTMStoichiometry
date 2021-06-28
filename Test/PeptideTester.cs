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
        private static readonly object[] _getModificationsTest =
        {
            new object[] {
                "KTEM[Common Variable:Oxidation on M]VSSVPAE[Metal:FeIII on E]NKSVLNEHQETSK",
                 new List<string>()
                 {
                     "Common Variable:Oxidation on M",
                     "Metal:FeIII on E"
                 },
                 new List<bool>()
                 {
                     true,
                     true
                 },
                 new List<string>()
                 {
                     "KTEMCommon Variable:Oxidation on MVSSVPAENKSVLNEHQETSK",
                     "KTEMVSSVPAEMetal:FeIII on ENKSVLNEHQETSK"
                 },
            },
            new object[] {
                "[Common Biological:Acetylation on X]RAEEPC[Common Fixed:Carbamidomethyl on C]APGAPSALGAQR",
                new List<string>()
                 {
                     "Common Biological:Acetylation on X",
                     "Common Fixed:Carbamidomethyl on C"
                 },
                 new List<bool>()
                 {
                     false,
                     true
                 },
                 new List<string>()
                 {
                     "Common Biological:Acetylation on XRAEEPCAPGAPSALGAQR",
                     "RAEEPCCommon Fixed:Carbamidomethyl on CAPGAPSALGAQR"
                 },
            },
        };

        [Test]
        //KTEM[Common Variable:Oxidation on M]VSSVPAE[Metal:Fe[III] on E]NKSVLNEHQETSK
        //[Common Biological:Acetylation on X]RAEEPC[Common Fixed:Carbamidomethyl on C]APGAPSALGAQR
        [TestCaseSource("_getModificationsTest")]
        public void Peptide_GetModifications_Pass(string seq, List<string> mods, List<bool> localized, List<string> localizedMods)
        {
            List<PostTranslationalModification> postTranslationalModifications = Peptide.GetModifications(seq);

            Assert.AreEqual(mods, postTranslationalModifications.Select(p => p.Modification));
            Assert.AreEqual(localized, postTranslationalModifications.Select(p => p.Localized));
            Assert.AreEqual(localizedMods, postTranslationalModifications.Select(p => p.ModificationInPeptideSequence));
        }

        //test that Peptide is reading in information from FlashLFQ lines and setting isUnique correctly
        [Test]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 0,
            "AGLGTGAAGGIGAGRTRAPSLASSSGSDK", "AGLGTGAAGGIGAGRTRAPSLASSSGSDK", "A0A087WPF7", "Auts2", "Mus musculus", 36, 0,
            "1DLC122419QE_ZD_Kidney_global_#4-calib", "Day 0", 0, true, false)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 1,
            "TPPSFPTPPPWLKPGELER", "TPPSFPTPPPWLKPGELER", "A0A087WPF7", "Auts2", "Mus musculus", 36, 0,
            "1DLC122419QE_ZD_Kidney_global_#4-calib", "Day 0", 0, true, false)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 2,
            "TPPTAALS[Common Biological:Phosphorylation on S]APPPLIS[Common Biological:Phosphorylation on S]TLGGR", "TPPTAALSAPPPLISTLGGR",
            "A0A087WPF7", "Auts2", "Mus musculus", 36, 0, "1DLC122419QE_ZD_Kidney_global_#4-calib", "Day 0", 0, true, false)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 2,
            "TPPTAALS[Common Biological:Phosphorylation on S]APPPLIS[Common Biological:Phosphorylation on S]TLGGR", "TPPTAALSAPPPLISTLGGR",
            "A0A087WPF7", "Auts2", "Mus musculus", 36, 6, "1DLC122419QE_ZD_Kidney_global_#15-calib", "Day 7", 8975790.34768643,  true, false)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 4,
            "ASVPSSEAGAWEVAASDIEPESRDRR", "ASVPSSEAGAWEVAASDIEPESRDRR", "A0A140LHF2", "Vsig10l2", "Mus musculus", 36, 35,
            "1DLC122919QE_ZD_Kidney_phos_#14-calib", "Day 7", 0, true, false)]

        public void Peptide_FlashLFQGetter_Pass(string filePath, string groupPath, int reqNumPepMeasurements, int lineToCompare, string sequence, string baseSequence,
            string proteinGroup, string geneName, string organism, int numberIntensities, int intensityColToCompare, string fileName,
            string groupID, double intensityVal, bool isUnique, bool detectedMinNum)
        {
            //check that valid entries in constuctor lead to object in getter
            List<Peptide> pepsInFile = PeptideReader.ReadTsv(filePath, groupPath, reqNumPepMeasurements, 5, "FlashLFQ");
            Assert.AreEqual(sequence, pepsInFile[lineToCompare].Sequence);
            Assert.AreEqual(baseSequence, pepsInFile[lineToCompare].BaseSeq);
            Assert.AreEqual(proteinGroup, pepsInFile[lineToCompare].ProteinGroup);
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

        //change only group file - everything else should match, so can check that all intensities are being read in correctly
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\twoGroupsPhosphoStudy.txt", 3)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPosphoSmallTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\twoGroupsPhosphoStudy.txt", 3)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPosphoLargeTest.tsv",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\twoGroupsPhosphoStudy.txt", 3)]
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

        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 0, 5)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 1, 2)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 2, 1)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 0)]
        public void Peptide_DetectedMinNum_Pass(string filePath, string groupPath, int reqNumPepMeasurements, int numPeptidesMeetMeasurementReq)
        {
            List<Peptide> pepsInFile = PeptideReader.ReadTsv(filePath, groupPath, reqNumPepMeasurements, 5, "FlashLFQ");
            Assert.AreEqual(numPeptidesMeetMeasurementReq, pepsInFile.Where(p => p.DetectedMinNum).Count());
        }

        private static readonly object[] _isUniqueLists =
        {
            new object[] { new List<Peptide> {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism", new List<Intensity>() { new Intensity("file", "group", 500) }, new List<string>() { "group" }, 0),
                new Peptide("Seq2", "Seq1", "Prot", "Gene", "Organism", new List<Intensity>() { new Intensity("file", "group", 500) }, new List<string>() { "group" }, 0),
                new Peptide("Seq3", "Seq3", "Prot", "Gene", "Organism", new List<Intensity>() { new Intensity("file", "group", 500) }, new List<string>() { "group" }, 0),
                new Peptide("Seq4", "Seq4", "Prot", "Gene", "Organism", new List<Intensity>() { new Intensity("file", "group", 500) }, new List<string>() { "group" }, 0),
                new Peptide("Seq5", "Seq5", "Prot", "Gene", "Organism", new List<Intensity>() { new Intensity("file", "group", 500) }, new List<string>() { "group" }, 0),
                new Peptide("Seq6", "Seq6", "Prot", "Gene", "Organism", new List<Intensity>() { new Intensity("file", "group", 500) }, new List<string>() { "group" }, 0)}, 
                new List<bool> {true, true, true, true, true, true}
            },
            new object[] { new List<Peptide> {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism", new List<Intensity>() { new Intensity("file", "group", 500) }, new List<string>() { "group" }, 0),
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism", new List<Intensity>() { new Intensity("file", "group", 500) }, new List<string>() { "group" }, 0),
                new Peptide("Seq1a", "Seq1", "Prot", "Gene", "Organism", new List<Intensity>() { new Intensity("file", "group", 500) }, new List<string>() { "group" }, 0),
                new Peptide("Seq2", "Seq2", "Prot", "Gene", "Organism", new List<Intensity>() { new Intensity("file", "group", 500) }, new List<string>() { "group" }, 0),
                new Peptide("Seq3", "Seq3", "Prot", "Gene", "Organism", new List<Intensity>() { new Intensity("file", "group", 500) }, new List<string>() { "group" }, 0),
                new Peptide("Seq2", "Seq4", "Prot", "Gene", "Organism", new List<Intensity>() { new Intensity("file", "group", 500) }, new List<string>() { "group" }, 0),
                new Peptide("Seq5", "Seq1", "Prot", "Gene", "Organism", new List<Intensity>() { new Intensity("file", "group", 500) }, new List<string>() { "group" }, 0),
                new Peptide("Seq6", "Seq6", "Prot", "Gene", "Organism", new List<Intensity>() { new Intensity("file", "group", 500) }, new List<string>() { "group" }, 0)},
                new List<bool> {false, false, true, false, true, false, true, true}
            }
        };

        /*

        [Test]

        [TestCaseSource("_isUniqueLists")]
        public void Peptide_IsUnique_Pass(List<Peptide> peps, List<bool> isUnique)
        {
            foreach (Peptide pep in peps)
            {
                pep.setIsUnique(peps);
            }
            Assert.AreEqual(isUnique, peps.Select(peps => peps.IsUnique));
        }
        */
    }
}
