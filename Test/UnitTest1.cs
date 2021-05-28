using NUnit.Framework;
using System;
using System.Collections.Generic;
using PTMStoichiometry20210414a;
using System.Linq;

namespace PTMStoichiometryTester20200415a
{
    //verification of program function
    [TestFixture]
    public class Tests
    {

        //Tests:

        //Clone: make sure that two identical imputs get identical results (maybe also import sorted and shuffled data -> should get equal output)

        //Mini data set: run through R/exccel: validate statictics

        //validate gets etc.
        //check pep1, pep2 NOT pep2, pep1 (maybe check that get same value if do calc both ways though) -> same for groups
        //check max limits (e.g. do not have more stoichiometries than samples in group * samples in group, or more peptides than started with, only using values with at least three stoichiometries)
        //validate all proteins being used are "useful"
        //validate all stoichiometires being used are useful
        //benjamini-Hochberg validation

        //want to have multiple test sets -> work through testing all functions (right now have phospho)

        //TODO: input & output validation
        

        //tiny test set (phospho): 5 peptides, 3 groups
        public static string filePathAllQuantPeptidesTinyTest = @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt";
        

        //small test set (phospho): 31 peptides, 3 groups
        public static string filePathAllQuantPeptidesSmallTest = @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPosphoSmallTest.txt";

        //large test set (phospho): 31 peptides, 3 groups
        public static string filePathAllQuantPeptidesLargeTest = @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPosphoLargeTest.tsv";

        //groups for Phospho data
        public static string filePathGroupsPhosphoTest = @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt";

        //test that PeptideReader reads in the correct number of Peptides
        [Test]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPosphoSmallTest.txt", 
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 31)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPosphoLargeTest.tsv",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 39869)]
        public void PeptideReaderTester(string filePath, string groupPath, int reqNumPepMeasurements, int peptideCount)
        {
            //check throws error if file references no good
            //var noFilePeptideReaderTest = PeptideReader.ReadTsv("", filePathGroupsOnePeptideTest);
            //Assert.Throws(, onePeptideReaderTest.Count);

            //small test: 31 peptides
            //List<Peptide> smallPeptideReaderTest = PeptideReader.ReadTsv(filePathAllQuantPeptidesSmallTest, filePathGroupsPhosphoTest);
            //Assert.AreEqual(31, smallPeptideReaderTest.Count);

            //large test: 39869 peptides
            //List<Peptide> largePeptideReaderTest = PeptideReader.ReadTsv(filePathAllQuantPeptidesLargeTest, filePathGroupsPhosphoTest);
            //Assert.AreEqual(39869, largePeptideReaderTest.Count);

            //tests that the correct number of peptides are being read in
            List<Peptide> pepsInFile = PeptideReader.ReadTsv(filePath, groupPath, reqNumPepMeasurements);
            Assert.AreEqual(peptideCount, pepsInFile.Count());

            //test group dictonaries & lists!
        }

        //test that Intensity is getting information correctly
        [Test]
        [TestCase("file.txt", "Control", 284.6, DetectionMS.MSMS)]
        [TestCase("file3.txt", "group 45", 3579.98456, DetectionMS.MS)]
        [TestCase("file3.raw", "group-42%AA", 0, DetectionMS.MSMSIdentifiedButNotQuantified)]
        [TestCase("file 4", "98%HA", 0, DetectionMS.NotDetected)]
        [TestCase("file_46_4.raw", "Disease 01", 0, DetectionMS.MSMSAmbiguousPeakfinding)]
        public void IntensityTester(string file, string group, double intensity, DetectionMS detection)
        {
            //check error thrown if entries invalid?

            //check valid entries in constuctor lead to object in getter
            Intensity TestIntensity = new Intensity(file, group, intensity, detection);
            Assert.AreEqual(file, TestIntensity.FileName);
            Assert.AreEqual(group, TestIntensity.GroupID);
            Assert.AreEqual(intensity, TestIntensity.IntensityVal);
            Assert.AreEqual(detection, TestIntensity.Detection);
        }

        //Test that stoichiometry is performing calculation correctly 
        //TODO: add pep:pep method, test that selects correct method
        [Test]
        [TestCase("file.txt", "Control", 284.6, DetectionMS.MSMS, 45737.7547, true)]
        [TestCase("file.txt", "Control", 2654754284.5686586, DetectionMS.MSMS, 45737.7547, true)]
        [TestCase("file3.txt", "group 45", 3579.98456, DetectionMS.MS, 27346.3456, true)]
        [TestCase("file3.txt", "group 45", 9843.43, DetectionMS.MS, 273443536.3456, true)]
        [TestCase("file3.raw", "group-42%AA", 0, DetectionMS.MSMSIdentifiedButNotQuantified, 3857.48, false)]
        [TestCase("file 4", "98%HA", 0, DetectionMS.NotDetected, 476.64, false)]
        [TestCase("file_46_4.raw", "Disease 01", 0, DetectionMS.MSMSAmbiguousPeakfinding, 1461636.42, false)]
        public void StoichiometryTester(string file, string group, double intensity, DetectionMS detection, double baseline, Boolean useful)
        {
            //check error thrown if entries invalid? (divide by 0 - check baseline > 0)

            //check valid entries in constuctor lead to object in getter
            Intensity Testintensity = new Intensity(file, group, intensity, detection);
            Stoichiometry TestStoichiometry = new Stoichiometry(Testintensity, baseline);
            Assert.AreEqual(Testintensity.IntensityVal/baseline, TestStoichiometry.StoichiometryVal);
            Assert.AreEqual(useful, TestStoichiometry.usefulStoichiometry);
        }

        //test that Peptide is reading in information from lines correctly and setting isUnique correctly
        [Test]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 0, 
            "AGLGTGAAGGIGAGRTRAPSLASSSGSDK", "AGLGTGAAGGIGAGRTRAPSLASSSGSDK", "A0A087WPF7", "Auts2", "Mus musculus", 36, 0,
            "1DLC122419QE_ZD_Kidney_global_#4-calib", "Day 0", 0, DetectionMS.NotDetected, true, false)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 1,
            "TPPSFPTPPPWLKPGELER", "TPPSFPTPPPWLKPGELER", "A0A087WPF7", "Auts2", "Mus musculus", 36, 0,
            "1DLC122419QE_ZD_Kidney_global_#4-calib", "Day 0", 0, DetectionMS.MSMSIdentifiedButNotQuantified, true, false)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 2,
            "TPPTAALS[Common Biological:Phosphorylation on S]APPPLIS[Common Biological:Phosphorylation on S]TLGGR", "TPPTAALSAPPPLISTLGGR", 
            "A0A087WPF7", "Auts2", "Mus musculus", 36, 0, "1DLC122419QE_ZD_Kidney_global_#4-calib", "Day 0", 0, DetectionMS.MSMSAmbiguousPeakfinding, true, false)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 2,
            "TPPTAALS[Common Biological:Phosphorylation on S]APPPLIS[Common Biological:Phosphorylation on S]TLGGR", "TPPTAALSAPPPLISTLGGR",
            "A0A087WPF7", "Auts2", "Mus musculus", 36, 6, "1DLC122419QE_ZD_Kidney_global_#15-calib", "Day 7", 8975790.34768643, DetectionMS.MS, true, false)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3, 4,
            "ASVPSSEAGAWEVAASDIEPESRDRR", "ASVPSSEAGAWEVAASDIEPESRDRR", "A0A140LHF2", "Vsig10l2", "Mus musculus", 36, 35,
            "1DLC122919QE_ZD_Kidney_phos_#14-calib", "Day 7", 0, DetectionMS.NotDetected, true, false)]

        public void PeptideTester(string filePath, string groupPath, int reqNumPepMeasurements, int lineToCompare, string sequence, string baseSequence, 
            string proteinGroup, string geneName, string organism, int numberIntensities, int intensityColToCompare, string fileName,
            string groupID, double intensityVal, DetectionMS detection, Boolean isUnique, Boolean detectedMinNum)
        {
            //check error thrown if entries invalid? (divide by 0 - check baseline > 0)

            //check valid entries in constuctor lead to object in getter
            List<Peptide> pepsInFile = PeptideReader.ReadTsv(filePath, groupPath, reqNumPepMeasurements);
            Assert.AreEqual(sequence, pepsInFile[lineToCompare].Sequence);
            Assert.AreEqual(baseSequence, pepsInFile[lineToCompare].BaseSeq);
            Assert.AreEqual(proteinGroup, pepsInFile[lineToCompare].ProteinGroup);
            Assert.AreEqual(geneName, pepsInFile[lineToCompare].GeneName);
            Assert.AreEqual(organism, pepsInFile[lineToCompare].Organism);
            Assert.AreEqual(numberIntensities, pepsInFile[lineToCompare].Intensities.Count);
            Assert.AreEqual(detectedMinNum, pepsInFile[lineToCompare].DetectedMinNum);

            //check setting isUnique
            Extensions.IncludeSharedPeptides(pepsInFile, true);
            Assert.AreEqual(isUnique, pepsInFile[lineToCompare].IsUnique);

            //check intensities
            Assert.AreEqual(fileName, pepsInFile[lineToCompare].Intensities[intensityColToCompare].FileName);
            Assert.AreEqual(groupID, pepsInFile[lineToCompare].Intensities[intensityColToCompare].GroupID);
            Assert.That(intensityVal, Is.EqualTo(pepsInFile[lineToCompare].Intensities[intensityColToCompare].IntensityVal).Within(0.001));
            Assert.AreEqual(detection, pepsInFile[lineToCompare].Intensities[intensityColToCompare].Detection);
        }

        //test that all shared peptides are being identified by IncludeSharedPeptidesTester
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinySharedPeptideTest.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 3,
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt", 2)]
        [Test]
        public void IncludeSharedPeptidesTester(string filePathSharedPeptides, string groupPath, int reqNumPepMeasurements, string filePathUniquePeptides, int numRazorPep)
        {
            List<Peptide> sharedPeps = PeptideReader.ReadTsv(filePathSharedPeptides, groupPath, reqNumPepMeasurements);
            List<Peptide> uniquePeps = PeptideReader.ReadTsv(filePathUniquePeptides, groupPath, reqNumPepMeasurements);

            //include razor peptides
            sharedPeps = Extensions.IncludeSharedPeptides(sharedPeps, true);
            Assert.AreNotEqual(uniquePeps.Count(), sharedPeps.Count());
            //remove razor peptides
            sharedPeps = Extensions.IncludeSharedPeptides(sharedPeps, false);
     
            Assert.AreEqual(uniquePeps.Count(), sharedPeps.Count());
        }
        


        [Test]
        public void BenjaminiHochbergTester(List<PairwiseCompairison> pairwiseComparisons, double alpha, double correctPVal)
        {
            Extensions.BenjaminiHochberg(pairwiseComparisons, alpha);
            //check that the p-values being set are correct
            //Assert.AreEqual(correctPVal, BHpval);
        }

        
        
        [Test]
        public void PairwiseCompairisonTester(Peptide Pep, List<Intensity> BaselinePepsIntensity, string G1, string G2, int minNumStoichiometries, 
            int numStoichiometriesGroupOne, int numStoichiometriesGroupTwo, int indexStoicToCheckG1, double stoichG1, int indexStoicToCheckG2, double stoichG2)
        {
            //check baseline case
            PairwiseCompairison TestPairwiseComparison = new PairwiseCompairison(Pep, BaselinePepsIntensity, G1, G2, minNumStoichiometries);
            Assert.AreEqual(Pep, TestPairwiseComparison.PeptideOne);
            Assert.AreEqual(G1, TestPairwiseComparison.GroupOne);
            Assert.AreEqual(G2, TestPairwiseComparison.GroupTwo);
            Assert.AreEqual(numStoichiometriesGroupOne, TestPairwiseComparison.PeptideStoichiometriesGroupOne.Count());
            Assert.AreEqual(numStoichiometriesGroupTwo, TestPairwiseComparison.PeptideStoichiometriesGroupTwo.Count());

            //check one stoichiometry in each group
            Assert.That(stoichG1, Is.EqualTo(TestPairwiseComparison.PeptideStoichiometriesGroupOne[indexStoicToCheckG1].StoichiometryVal).Within(0.001));
            Assert.That(stoichG2, Is.EqualTo(TestPairwiseComparison.PeptideStoichiometriesGroupTwo[indexStoicToCheckG2].StoichiometryVal).Within(0.001));
            //check 0 and NA values are being removed 
            Assert.IsEmpty(TestPairwiseComparison.PeptideStoichiometriesGroupOne.Where(p => !p.usefulStoichiometry));
            Assert.IsEmpty(TestPairwiseComparison.PeptideStoichiometriesGroupTwo.Where(p => !p.usefulStoichiometry));

            //check peptide:peptide case
        }







































        
        //test baseline peptides - make case where clearly covary & where clearly don't covary
        //through test case data sets check proteingroup, pairwisecompairison

        //have a tester for extensions methods: BH, setting unique peptides

        //tester for writer?


        //baseline peptides tester (this will test covaraince as well)


        [Test]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\AllQuantifiedPeptidesPhosphoTinyTest.txt", 
            @"C:\Users\KAP\source\repos\PTMStoichiometryTester20200415a\TestData\groupsPhosphoStudy.txt", 0, 1, "TPPSFPTPPPWLKPGELER", "A0A087WPF7",
            "InsufficientPeptides", 3, null, null, null, null)]

        public void ProteinGroupTester(string filePath, string groupPath, int reqNumPepMeasurements, int reqNumUnmodPeptides, int reqNumModPeptides, int reqNumOfPepeptides,
        Boolean useBaselinePeptides, int reqNumBaselinePeptides, double correlationCutOff, Boolean compareUnmod, int minNumStoichiometries, 
        int lineToCompare, int IndexOfPeptideToCheck, string PeptideSeq, string ProteinName, string useProt, int NumPeptidesInProtein, 
        int IndexOfBaselinePeptideToCheck, string BaselinePeptideBaseSeq, int NumBaselinePeptides, int NumPairwiseCompairisons)
        {

            Dictionary<string, string> groups = PeptideReader.GetGroups(groupPath);
            List<string> groupsList = PeptideReader.GetGroupList(groupPath);
            List<Peptide> testPeptide = PeptideReader.ReadTsv(filePath, groupPath, reqNumPepMeasurements);
            Extensions.IncludeSharedPeptides(testPeptide, false);

            //group peptides by protein
            var proteins = testPeptide.Select(p => p.ProteinGroup).Distinct().ToArray();
            List<ProteinGroup> testProt = new List<ProteinGroup>();
            for (int i = 0; i < proteins.Length; i++)
            {
                testProt.Add(new ProteinGroup(proteins[i], testPeptide, groupsList, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides,
                    useBaselinePeptides, reqNumBaselinePeptides, correlationCutOff, compareUnmod, minNumStoichiometries));
            }
            testProt.OrderBy(p => p.ProteinName);
            
            Assert.AreEqual(ProteinName, testProt[lineToCompare].ProteinName);
            Assert.AreEqual(NumPeptidesInProtein, testProt[lineToCompare].PeptidesInProtein.Count());
            Assert.AreEqual(NumPeptidesInProtein, testProt[lineToCompare].NumPeptidesInProtein);
            Assert.AreEqual(PeptideSeq, testProt[lineToCompare].PeptidesInProtein[IndexOfPeptideToCheck].BaseSeq);
            Assert.AreEqual(useProt, testProt[lineToCompare].useProt);
            
            if (testProt[lineToCompare].useProt)
            {
                Assert.AreEqual(NumBaselinePeptides, testProt[lineToCompare].BaselinePeptides.Count());
                Assert.AreEqual(BaselinePeptideBaseSeq, testProt[lineToCompare].BaselinePeptides[IndexOfBaselinePeptideToCheck].BaseSeq);
                if (testProt[lineToCompare].BaselinePeptides.Count() > 2)
                {
                    Assert.AreEqual(NumPairwiseCompairisons, testProt[lineToCompare].ProteinPairwiseComparisons.Count());
                }
                else
                {
                    Assert.IsNull(testProt[lineToCompare].ProteinPairwiseComparisons);
                }
            }
            else
            {
                Assert.IsNull(testProt[lineToCompare].BaselinePeptides);
            }

            





            //check different ways of narrowing down the protein groups
            testProt = testProt.Where(p => p.ProteinPairwiseComparisons != null).ToList();
            Extensions.CalcCorrectedPValue(testProt, false);
            List<ProteinGroup> ProteinsToUse = testProt.Where(p => p.useProt).ToList();
           
        }





       /*
        [Test]
        public void TinyPhosphoTestCase()
        {
            //run of program on small phospho data set: 31 peptides, 3 groups
            //read in peptide data
            List<string> groupsList = PeptideReader.GetGroupList(filePathGroupsPhosphoTest);
            List<Peptide> tinyPeptideTest = PeptideReader.ReadTsv(filePathAllQuantPeptidesTinyTest, filePathGroupsPhosphoTest);

            //group peptides by protein
            var proteins = tinyPeptideTest.Select(p => p.ProteinGroup).Distinct().ToArray();
            List<ProteinGroup> testProt = new List<ProteinGroup>();
            for (int i = 0; i < proteins.Length; i++)
            {
                testProt.Add(new ProteinGroup(proteins[i], tinyPeptideTest, groupsList));
            }
            List<ProteinGroup> proteinsToUse = testProt.Where(p => p.useProt).ToList();

            //validate Peptide data input correctly: check correct length and all Peptides
            Assert.AreEqual(5, tinyPeptideTest.Count());
            Assert.AreEqual("AGLGTGAAGGIGAGRTRAPSLASSSGSDK", tinyPeptideTest[0].Sequence);
            Assert.AreEqual("AGLGTGAAGGIGAGRTRAPSLASSSGSDK", tinyPeptideTest[0].BaseSeq);
            Assert.AreEqual("A0A087WPF7", tinyPeptideTest[0].ProteinGroup);
            Assert.AreEqual("Auts2", tinyPeptideTest[0].GeneName);
            Assert.AreEqual("Mus musculus", tinyPeptideTest[0].Organism);
            Assert.AreEqual(36, tinyPeptideTest[0].Intensities.Count);

            Assert.AreEqual("TPPSFPTPPPWLKPGELER", tinyPeptideTest[1].Sequence);
            Assert.AreEqual("TPPSFPTPPPWLKPGELER", tinyPeptideTest[1].BaseSeq);
            Assert.AreEqual("A0A087WPF7", tinyPeptideTest[1].ProteinGroup);
            Assert.AreEqual("Auts2", tinyPeptideTest[1].GeneName);
            Assert.AreEqual("Mus musculus", tinyPeptideTest[1].Organism);
            Assert.AreEqual(36, tinyPeptideTest[1].Intensities.Count);

            Assert.AreEqual("TPPTAALS[Common Biological:Phosphorylation on S]APPPLIS[Common Biological:Phosphorylation on S]TLGGR", tinyPeptideTest[2].Sequence);
            Assert.AreEqual("TPPTAALSAPPPLISTLGGR", tinyPeptideTest[2].BaseSeq);
            Assert.AreEqual("A0A087WPF7", tinyPeptideTest[2].ProteinGroup);
            Assert.AreEqual("Auts2", tinyPeptideTest[2].GeneName);
            Assert.AreEqual("Mus musculus", tinyPeptideTest[2].Organism);
            Assert.AreEqual(36, tinyPeptideTest[2].Intensities.Count);

            Assert.AreEqual("PSLDQIYTQFK", tinyPeptideTest[3].Sequence);
            Assert.AreEqual("PSLDQIYTQFK", tinyPeptideTest[3].BaseSeq);
            Assert.AreEqual("A0A0U1RPR8", tinyPeptideTest[3].ProteinGroup);
            Assert.AreEqual("Gucy2d", tinyPeptideTest[3].GeneName);
            Assert.AreEqual("Mus musculus", tinyPeptideTest[3].Organism);
            Assert.AreEqual(36, tinyPeptideTest[3].Intensities.Count);

            Assert.AreEqual("ASVPSSEAGAWEVAASDIEPESRDRR", tinyPeptideTest[4].Sequence);
            Assert.AreEqual("ASVPSSEAGAWEVAASDIEPESRDRR", tinyPeptideTest[4].BaseSeq);
            Assert.AreEqual("A0A140LHF2", tinyPeptideTest[4].ProteinGroup);
            Assert.AreEqual("Vsig10l2", tinyPeptideTest[4].GeneName);
            Assert.AreEqual("Mus musculus", tinyPeptideTest[4].Organism);
            Assert.AreEqual(36, tinyPeptideTest[4].Intensities.Count);

            //Validate intensity values: check length (checked above), first Intensity, last Intensity and select middle Intensities 
            Assert.AreEqual("1DLC122419QE_ZD_Kidney_global_#4-calib", tinyPeptideTest[0].Intensities[0].FileName);
            Assert.AreEqual("Day 0", tinyPeptideTest[0].Intensities[0].GroupID);
            Assert.AreEqual(0, tinyPeptideTest[0].Intensities[0].IntensityVal);
            Assert.AreEqual(DetectionMS.NotDetected, tinyPeptideTest[0].Intensities[0].Detection);

            Assert.AreEqual("1DLC122419QE_ZD_Kidney_global_#4-calib", tinyPeptideTest[1].Intensities[0].FileName);
            Assert.AreEqual("Day 0", tinyPeptideTest[1].Intensities[0].GroupID);
            Assert.AreEqual(0, tinyPeptideTest[1].Intensities[0].IntensityVal);
            Assert.AreEqual(DetectionMS.MSMSIdentifiedButNotQuantified, tinyPeptideTest[1].Intensities[0].Detection);

            Assert.AreEqual("1DLC122419QE_ZD_Kidney_global_#4-calib", tinyPeptideTest[2].Intensities[0].FileName);
            Assert.AreEqual("Day 0", tinyPeptideTest[2].Intensities[0].GroupID);
            Assert.AreEqual(0, tinyPeptideTest[2].Intensities[0].IntensityVal);
            Assert.AreEqual(DetectionMS.MSMSAmbiguousPeakfinding, tinyPeptideTest[2].Intensities[0].Detection);

            Assert.AreEqual("1DLC122419QE_ZD_Kidney_global_#10-calib", tinyPeptideTest[0].Intensities[3].FileName);
            Assert.AreEqual("Day 2", tinyPeptideTest[0].Intensities[3].GroupID);
            Assert.That(2201586.51483236, Is.EqualTo(tinyPeptideTest[0].Intensities[3].IntensityVal).Within(0.001));
            Assert.AreEqual(DetectionMS.MSMS, tinyPeptideTest[0].Intensities[3].Detection);

            Assert.AreEqual("1DLC122419QE_ZD_Kidney_global_#15-calib", tinyPeptideTest[2].Intensities[6].FileName);
            Assert.AreEqual("Day 7", tinyPeptideTest[2].Intensities[6].GroupID);
            Assert.That(8975790.34768643, Is.EqualTo(tinyPeptideTest[2].Intensities[6].IntensityVal).Within(0.001));
            Assert.AreEqual(DetectionMS.MS, tinyPeptideTest[2].Intensities[6].Detection);

            Assert.AreEqual("1DLC122919QE_ZD_Kidney_phos_#14-calib", tinyPeptideTest[4].Intensities[35].FileName);
            Assert.AreEqual("Day 7", tinyPeptideTest[4].Intensities[35].GroupID);
            Assert.AreEqual(0, tinyPeptideTest[4].Intensities[35].IntensityVal);
            Assert.AreEqual(DetectionMS.NotDetected, tinyPeptideTest[4].Intensities[35].Detection);
            
            //validate protein group significance: p-values calulcated in.... 
            //check correct number of proteins (only one protein with mulitple peptides in set)
            Assert.AreEqual(1, proteinsToUse.Count());
            
            //verify the ProteinGroup
            Assert.AreEqual("A0A087WPF7", proteinsToUse[0].ProteinName);
            Assert.AreEqual("AGLGTGAAGGIGAGRTRAPSLASSSGSDK", proteinsToUse[0].PeptidesInProtein[0].Sequence);
            Assert.AreEqual("TPPSFPTPPPWLKPGELER", proteinsToUse[0].PeptidesInProtein[1].Sequence);
            Assert.AreEqual("TPPTAALS[Common Biological:Phosphorylation on S]APPPLIS[Common Biological:Phosphorylation on S]TLGGR", proteinsToUse[0].PeptidesInProtein[2].Sequence);
            Assert.AreEqual(3, proteinsToUse[0].NumPeptidesInProtein);
            Assert.AreEqual("modandunmod", proteinsToUse[0].useProt);
            Assert.AreEqual(1, proteinsToUse[0].ProteinPairwiseComparisons.Count()); // pep1/pep2, pep1/pep3, pep2/pep3 each 3x (bc there are 3 ways we can pair groups), -> filter for both dist at least 3 stoic -> only one
            //TODO: check BH p-value
            //Assert.That(0, Is.EqualTo(proteinsToUse[0].CorrectedMWPValue).Within(0.001));

            //verify the pairwise comparisons: check correct peptide & group pairs (no duplicates or self pairs)
            //AGLGTGAAGGIGAGRTRAPSLASSSGSDK, TPPTAALS[Common Biological:Phosphorylation on S]APPPLIS[Common Biological:Phosphorylation on S]TLGGR, Day 2, Day 7
            //Assert.AreEqual(tinyPeptideTest[0], proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideOne);
            //Assert.AreEqual(tinyPeptideTest[2], proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideTwo);
            Assert.AreEqual(groupsList[0], proteinsToUse[0].ProteinPairwiseComparisons[0].GroupOne);
            Assert.AreEqual(groupsList[1], proteinsToUse[0].ProteinPairwiseComparisons[0].GroupTwo);
            Assert.AreEqual(6, proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupOne.Count());
            Assert.AreEqual(8, proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupTwo.Count());
            Assert.AreEqual(14, proteinsToUse[0].ProteinPairwiseComparisons[0].MWStat);
            Assert.That(0.2284, Is.EqualTo(proteinsToUse[0].ProteinPairwiseComparisons[0].MWPVal).Within(0.05));
            */
            /*
            //verify the pairwise comparisons: check correct peptide & group pairs (no duplicates or self pairs)
            //AGLGTGAAGGIGAGRTRAPSLASSSGSDK, TPPSFPTPPPWLKPGELER, Day 2, Day 7
            Assert.AreEqual(tinyPeptideTest[0], proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideOne);
            Assert.AreEqual(tinyPeptideTest[1], proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideTwo);
            Assert.AreEqual(groupsList[0], proteinsToUse[0].ProteinPairwiseComparisons[0].GroupOne);
            Assert.AreEqual(groupsList[1], proteinsToUse[0].ProteinPairwiseComparisons[0].GroupTwo);
            Assert.AreEqual(0, proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupOne.Count());
            Assert.AreEqual(4, proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupTwo.Count());

            //AGLGTGAAGGIGAGRTRAPSLASSSGSDK, TPPSFPTPPPWLKPGELER, Day 2, Day 0
            Assert.AreEqual(tinyPeptideTest[0], proteinsToUse[0].ProteinPairwiseComparisons[1].PeptideOne);
            Assert.AreEqual(tinyPeptideTest[1], proteinsToUse[0].ProteinPairwiseComparisons[1].PeptideTwo);
            Assert.AreEqual(groupsList[0], proteinsToUse[0].ProteinPairwiseComparisons[1].GroupOne);
            Assert.AreEqual(groupsList[2], proteinsToUse[0].ProteinPairwiseComparisons[1].GroupTwo);
            Assert.AreEqual(0, proteinsToUse[0].ProteinPairwiseComparisons[1].PeptideStoichiometriesGroupOne.Count());
            Assert.AreEqual(0, proteinsToUse[0].ProteinPairwiseComparisons[1].PeptideStoichiometriesGroupTwo.Count());

            //AGLGTGAAGGIGAGRTRAPSLASSSGSDK, TPPSFPTPPPWLKPGELER, Day 7, Day 0
            Assert.AreEqual(tinyPeptideTest[0], proteinsToUse[0].ProteinPairwiseComparisons[2].PeptideOne);
            Assert.AreEqual(tinyPeptideTest[1], proteinsToUse[0].ProteinPairwiseComparisons[2].PeptideTwo);
            Assert.AreEqual(groupsList[1], proteinsToUse[0].ProteinPairwiseComparisons[2].GroupOne);
            Assert.AreEqual(groupsList[2], proteinsToUse[0].ProteinPairwiseComparisons[2].GroupTwo);
            Assert.AreEqual(4, proteinsToUse[0].ProteinPairwiseComparisons[2].PeptideStoichiometriesGroupOne.Count());
            Assert.AreEqual(0, proteinsToUse[0].ProteinPairwiseComparisons[2].PeptideStoichiometriesGroupTwo.Count());

            //AGLGTGAAGGIGAGRTRAPSLASSSGSDK, TPPTAALS[Common Biological:Phosphorylation on S]APPPLIS[Common Biological:Phosphorylation on S]TLGGR, Day 2, Day 7
            Assert.AreEqual(tinyPeptideTest[0], proteinsToUse[0].ProteinPairwiseComparisons[3].PeptideOne);
            Assert.AreEqual(tinyPeptideTest[2], proteinsToUse[0].ProteinPairwiseComparisons[3].PeptideTwo);
            Assert.AreEqual(groupsList[0], proteinsToUse[0].ProteinPairwiseComparisons[3].GroupOne);
            Assert.AreEqual(groupsList[1], proteinsToUse[0].ProteinPairwiseComparisons[3].GroupTwo);
            Assert.AreEqual(6, proteinsToUse[0].ProteinPairwiseComparisons[3].PeptideStoichiometriesGroupOne.Count());
            Assert.AreEqual(8, proteinsToUse[0].ProteinPairwiseComparisons[3].PeptideStoichiometriesGroupTwo.Count());

            //AGLGTGAAGGIGAGRTRAPSLASSSGSDK, TPPTAALS[Common Biological:Phosphorylation on S]APPPLIS[Common Biological:Phosphorylation on S]TLGGR, Day 2, Day 0
            Assert.AreEqual(tinyPeptideTest[0], proteinsToUse[0].ProteinPairwiseComparisons[4].PeptideOne);
            Assert.AreEqual(tinyPeptideTest[2], proteinsToUse[0].ProteinPairwiseComparisons[4].PeptideTwo);
            Assert.AreEqual(groupsList[0], proteinsToUse[0].ProteinPairwiseComparisons[4].GroupOne);
            Assert.AreEqual(groupsList[2], proteinsToUse[0].ProteinPairwiseComparisons[4].GroupTwo);
            Assert.AreEqual(6, proteinsToUse[0].ProteinPairwiseComparisons[4].PeptideStoichiometriesGroupOne.Count());
            Assert.AreEqual(3, proteinsToUse[0].ProteinPairwiseComparisons[4].PeptideStoichiometriesGroupTwo.Count());

            //AGLGTGAAGGIGAGRTRAPSLASSSGSDK, TPPTAALS[Common Biological:Phosphorylation on S]APPPLIS[Common Biological:Phosphorylation on S]TLGGR, Day 7, Day 0
            Assert.AreEqual(tinyPeptideTest[0], proteinsToUse[0].ProteinPairwiseComparisons[5].PeptideOne);
            Assert.AreEqual(tinyPeptideTest[2], proteinsToUse[0].ProteinPairwiseComparisons[5].PeptideTwo);
            Assert.AreEqual(groupsList[1], proteinsToUse[0].ProteinPairwiseComparisons[5].GroupOne);
            Assert.AreEqual(groupsList[2], proteinsToUse[0].ProteinPairwiseComparisons[5].GroupTwo);
            Assert.AreEqual(8, proteinsToUse[0].ProteinPairwiseComparisons[5].PeptideStoichiometriesGroupOne.Count());
            Assert.AreEqual(3, proteinsToUse[0].ProteinPairwiseComparisons[5].PeptideStoichiometriesGroupTwo.Count());

            //TPPSFPTPPPWLKPGELER, TPPTAALS[Common Biological:Phosphorylation on S]APPPLIS[Common Biological:Phosphorylation on S]TLGGR, Day 2, Day 7
            Assert.AreEqual(tinyPeptideTest[1], proteinsToUse[0].ProteinPairwiseComparisons[6].PeptideOne);
            Assert.AreEqual(tinyPeptideTest[2], proteinsToUse[0].ProteinPairwiseComparisons[6].PeptideTwo);
            Assert.AreEqual(groupsList[0], proteinsToUse[0].ProteinPairwiseComparisons[6].GroupOne);
            Assert.AreEqual(groupsList[1], proteinsToUse[0].ProteinPairwiseComparisons[6].GroupTwo);
            Assert.AreEqual(0, proteinsToUse[0].ProteinPairwiseComparisons[6].PeptideStoichiometriesGroupOne.Count());
            Assert.AreEqual(2, proteinsToUse[0].ProteinPairwiseComparisons[6].PeptideStoichiometriesGroupTwo.Count());

            //TPPSFPTPPPWLKPGELER, TPPTAALS[Common Biological:Phosphorylation on S]APPPLIS[Common Biological:Phosphorylation on S]TLGGR, Day 2, Day 0
            Assert.AreEqual(tinyPeptideTest[1], proteinsToUse[0].ProteinPairwiseComparisons[7].PeptideOne);
            Assert.AreEqual(tinyPeptideTest[2], proteinsToUse[0].ProteinPairwiseComparisons[7].PeptideTwo);
            Assert.AreEqual(groupsList[0], proteinsToUse[0].ProteinPairwiseComparisons[7].GroupOne);
            Assert.AreEqual(groupsList[2], proteinsToUse[0].ProteinPairwiseComparisons[7].GroupTwo);
            Assert.AreEqual(0, proteinsToUse[0].ProteinPairwiseComparisons[7].PeptideStoichiometriesGroupOne.Count());
            Assert.AreEqual(0, proteinsToUse[0].ProteinPairwiseComparisons[7].PeptideStoichiometriesGroupTwo.Count());

            //TPPSFPTPPPWLKPGELER, TPPTAALS[Common Biological:Phosphorylation on S]APPPLIS[Common Biological:Phosphorylation on S]TLGGR, Day 7, Day 0
            Assert.AreEqual(tinyPeptideTest[1], proteinsToUse[0].ProteinPairwiseComparisons[8].PeptideOne);
            Assert.AreEqual(tinyPeptideTest[2], proteinsToUse[0].ProteinPairwiseComparisons[8].PeptideTwo);
            Assert.AreEqual(groupsList[1], proteinsToUse[0].ProteinPairwiseComparisons[8].GroupOne);
            Assert.AreEqual(groupsList[2], proteinsToUse[0].ProteinPairwiseComparisons[8].GroupTwo);
            Assert.AreEqual(2, proteinsToUse[0].ProteinPairwiseComparisons[8].PeptideStoichiometriesGroupOne.Count());
            Assert.AreEqual(0, proteinsToUse[0].ProteinPairwiseComparisons[8].PeptideStoichiometriesGroupTwo.Count());
            Assert.IsNull(proteinsToUse[0].ProteinPairwiseComparisons[8].MWStat);
            Assert.IsNull(proteinsToUse[0].ProteinPairwiseComparisons[8].MWPVal);
            */
            /*
            //verify the stoichiometries: checked lengths above - here checking values
            Assert.That(1.000212167, Is.EqualTo(proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupOne[0].StoichiometryVal).Within(0.0001));
            Assert.That(6.692823891, Is.EqualTo(proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupOne[1].StoichiometryVal).Within(0.0001));
            Assert.That(0.007443501, Is.EqualTo(proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupOne[2].StoichiometryVal).Within(0.0001));
            Assert.That(0.049807477, Is.EqualTo(proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupOne[3].StoichiometryVal).Within(0.0001));
            Assert.That(0.040703639, Is.EqualTo(proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupOne[4].StoichiometryVal).Within(0.0001));
            Assert.That(0.272364501, Is.EqualTo(proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupOne[5].StoichiometryVal).Within(0.0001));

            Assert.That(0.155590839, Is.EqualTo(proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupTwo[0].StoichiometryVal).Within(0.0001));
            Assert.That(58.23002238, Is.EqualTo(proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupTwo[1].StoichiometryVal).Within(0.0001));
            Assert.That(0.094474158, Is.EqualTo(proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupTwo[2].StoichiometryVal).Within(0.0001));
            Assert.That(35.35704522, Is.EqualTo(proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupTwo[3].StoichiometryVal).Within(0.0001));
            Assert.That(0.008434811, Is.EqualTo(proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupTwo[4].StoichiometryVal).Within(0.0001));
            Assert.That(3.15673613, Is.EqualTo(proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupTwo[5].StoichiometryVal).Within(0.0001));
            Assert.That(0.568020241, Is.EqualTo(proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupTwo[6].StoichiometryVal).Within(0.0001));
            Assert.That(212.5821259, Is.EqualTo(proteinsToUse[0].ProteinPairwiseComparisons[0].PeptideStoichiometriesGroupTwo[7].StoichiometryVal).Within(0.0001));
        }
            */
   

        /*
        [Test]
        public void SmallPhosphoTestCase()
        {
            //run of program on small phospho data set: 31 peptides, 3 groups
            //read in peptide data
            List<string> groupsList = PeptideReader.GetGroupList(filePathGroupsPhosphoTest);
            List<Peptide> testPeptide = PeptideReader.ReadTsv(filePathAllQuantPeptidesSmallTest, filePathGroupsPhosphoTest);

            //group peptides by protein
            var proteins = testPeptide.Select(p => p.ProteinGroup).Distinct().ToArray();
            List<ProteinGroup> testProt = new List<ProteinGroup>();
            for (int i = 0; i < proteins.Length; i++)
            {
                testProt.Add(new ProteinGroup(proteins[i], testPeptide, groupsList));
            }
            List<ProteinGroup> proteinsToUse = testProt.Where(p => p.useProt).ToList();

            //validate Peptide data input correctly

            //validate protein group significance: p-values calulcated in.... 

            //verify the pairwise comparisons: no pep1=pep2, no group1=group2, no pep1,pep2, then pep2,pep1...
        }
        */

        /*
        [Test]
        public void LargePhosphoTestCase()
        {
            //run of program on large phospho data set: 39869 peptides, 3 groups
            //read in peptide data
            List<string> groupsList = PeptideReader.GetGroupList(filePathGroupsPhosphoTest);
            List<Peptide> largePeptideTest = PeptideReader.ReadTsv(filePathAllQuantPeptidesLargeTest, filePathGroupsPhosphoTest);

            //group peptides by protein
            var proteins = largePeptideTest.Select(p => p.ProteinGroup).Distinct().ToArray();
            List<ProteinGroup> testProt = new List<ProteinGroup>();
            for (int i = 0; i < proteins.Length; i++)
            {
                testProt.Add(new ProteinGroup(proteins[i], largePeptideTest, groupsList));
            }
            List<ProteinGroup> proteinsToUse = testProt.Where(p => p.useProt).ToList();

            //validate Peptide data input correctly: check first, last Peptide, and select middle Peptides 
            Assert.AreEqual("AGLGTGAAGGIGAGRTRAPSLASSSGSDK", largePeptideTest[0].Sequence);
            Assert.AreEqual("AGLGTGAAGGIGAGRTRAPSLASSSGSDK", largePeptideTest[0].BaseSeq);
            Assert.AreEqual("A0A087WPF7", largePeptideTest[0].ProteinGroup);
            Assert.AreEqual("Auts2", largePeptideTest[0].GeneName);
            Assert.AreEqual("Mus musculus", largePeptideTest[0].Organism);
            Assert.AreEqual(36, largePeptideTest[0].Intensities.Count);

            Assert.AreEqual("KEELALLDQAAPM[Common Variable:Oxidation on M]ETNVTIKDLR", largePeptideTest[39868].Sequence);
            Assert.AreEqual("KEELALLDQAAPMETNVTIKDLR", largePeptideTest[39868].BaseSeq);
            Assert.AreEqual("W8DXL4", largePeptideTest[39868].ProteinGroup);
            Assert.AreEqual("Lrit3", largePeptideTest[39868].GeneName);
            Assert.AreEqual("Mus musculus", largePeptideTest[39868].Organism);
            Assert.AreEqual(36, largePeptideTest[39868].Intensities.Count);

            Assert.AreEqual("VSHALAEGLGVIAC[Common Fixed:Carbamidomethyl on C]IGEK", largePeptideTest[7094].Sequence);
            Assert.AreEqual("VSHALAEGLGVIACIGEK", largePeptideTest[7094].BaseSeq);
            Assert.AreEqual("P17751", largePeptideTest[7094].ProteinGroup);
            Assert.AreEqual("Tpi1", largePeptideTest[7094].GeneName);
            Assert.AreEqual("Mus musculus", largePeptideTest[7094].Organism);
            Assert.AreEqual(36, largePeptideTest[7094].Intensities.Count);

            Assert.AreEqual("VFLEELPPAT[UniProt:Phosphothreonine on T]PS[UniProt:Phosphoserine on S]PR", largePeptideTest[17994].Sequence);
            Assert.AreEqual("VFLEELPPATPSPR", largePeptideTest[17994].BaseSeq);
            Assert.AreEqual("Q60825", largePeptideTest[17994].ProteinGroup);
            Assert.AreEqual("Slc34a1", largePeptideTest[17994].GeneName);
            Assert.AreEqual("Mus musculus", largePeptideTest[17994].Organism);
            Assert.AreEqual(36, largePeptideTest[17994].Intensities.Count);

            //Validate intensity values: check first, last Peptide, and select middle Peptides 
            Assert.AreEqual("1DLC122419QE_ZD_Kidney_global_#4-calib", largePeptideTest[0].Intensities[0].FileName);
            Assert.AreEqual("Day 0", largePeptideTest[0].Intensities[0].GroupID);
            Assert.AreEqual(0, largePeptideTest[0].Intensities[0].IntensityVal);
            Assert.AreEqual(DetectionMS.NotDetected, largePeptideTest[0].Intensities[0].Detection);

            Assert.AreEqual("1DLC122919QE_ZD_Kidney_phos_#14-calib", largePeptideTest[211].Intensities[35].FileName);
            Assert.AreEqual("Day 7", largePeptideTest[211].Intensities[35].GroupID);
            Assert.AreEqual(0, largePeptideTest[211].Intensities[35].IntensityVal);
            Assert.AreEqual(DetectionMS.MSMSAmbiguousPeakfinding, largePeptideTest[211].Intensities[35].Detection);

            Assert.AreEqual("1DLC122419QE_ZD_Kidney_global_#9-calib", largePeptideTest[27].Intensities[4].FileName);
            Assert.AreEqual("Day 2", largePeptideTest[27].Intensities[4].GroupID);
            Assert.AreEqual(356346.1856, largePeptideTest[27].Intensities[4].IntensityVal);
            Assert.AreEqual(DetectionMS.MSMS, largePeptideTest[27].Intensities[4].Detection);

            Assert.AreEqual("1DLC122419QE_ZD_Kidney_global_#6-calib", largePeptideTest[27031].Intensities[11].FileName);
            Assert.AreEqual("Day 0", largePeptideTest[27031].Intensities[11].GroupID);
            Assert.AreEqual(0, largePeptideTest[27031].Intensities[11].IntensityVal);
            Assert.AreEqual(DetectionMS.MSMSIdentifiedButNotQuantified, largePeptideTest[27031].Intensities[11].Detection);

            //validate protein group significance: p-values calulcated in.... 
        }
        */


    }
}