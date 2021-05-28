using NUnit.Framework;
using System;
using System.Collections.Generic;
using PTMStoichiometry20210414a;
using System.Linq;

namespace PTMStoichiometryTester20200415a
{
    [TestFixture]
    class PairwiseCompairisonTester
    {
        
        private static readonly object[] _calcStoichiometryBaselineCase =
        {
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file", "group1", 0, DetectionMS.MSMSIdentifiedButNotQuantified),
                        new Intensity("file", "group2", 500, DetectionMS.MS) },
                    new List<string>() { "group1", "group2" }, 1),
                new List<Intensity> {
                new Intensity("file", "group1", 1000, DetectionMS.MS),
                new Intensity("file", "group2", 500, DetectionMS.MS)},
                "group1", "group2",
                3,
                new List<double> { },
                new List<double> { 1 }
            },
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism", 
                    new List<Intensity>() { 
                        new Intensity("file", "group1", 500, DetectionMS.MS), 
                        new Intensity("file", "group2", 500, DetectionMS.MS) }, 
                    new List<string>() { "group1", "group2" }, 1),
                new List<Intensity> {
                new Intensity("file", "group1", 1000, DetectionMS.MS),
                new Intensity("file", "group2", 250, DetectionMS.MS)},
                "group1", "group2",
                3,
                new List<double> { 0.5 },
                new List<double> { 2 }
            },
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism", 
                    new List<Intensity>() { 
                        new Intensity("file", "group1", 756019.37596012, DetectionMS.MS), 
                        new Intensity("file", "group2", 234667.023876, DetectionMS.MS) },
                    new List<string>() { "group1", "group2" }, 1),
                new List<Intensity> {
                new Intensity("file", "group1", 23974.23487, DetectionMS.MS),
                new Intensity("file", "group2", 807934.23874, DetectionMS.MS)},
                "group1", "group2",
                3,
                new List<double> { 756019.37596012/23974.23487 },
                new List<double> { 234667.023876/807934.23874 }
            },
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism", 
                    new List<Intensity>() { 
                        new Intensity("file", "group1", 756019.37596012, DetectionMS.MS), 
                        new Intensity("file", "group1", 234667.023876, DetectionMS.MS),
                        new Intensity("file", "group2", 38947.830453, DetectionMS.MS),
                        new Intensity("file", "group2", 8907235.0893645, DetectionMS.MS),
                        new Intensity("file", "group2", 0, DetectionMS.NotDetected),
                        new Intensity("file", "group3", 8907235.0893645, DetectionMS.MS),
                        new Intensity("file", "group3", 3409750.2309, DetectionMS.NotDetected)
                    },
                    new List<string>() { "group1", "group2" }, 1),
                new List<Intensity> {
                new Intensity("file", "group1", 23974.23487, DetectionMS.MS),
                new Intensity("file", "group1", 807934.23874, DetectionMS.MS),
                new Intensity("file", "group1", 23974.23487, DetectionMS.MS),
                new Intensity("file", "group2", 89754.37958, DetectionMS.MS),
                new Intensity("file", "group2", 357698.893275, DetectionMS.MS),
                new Intensity("file", "group2", 32587.4569384, DetectionMS.MS),
                new Intensity("file", "group3", 89754.37958, DetectionMS.MS),
                new Intensity("file", "group3", 357698.893275, DetectionMS.MS),
                new Intensity("file", "group3", 32587.4569384, DetectionMS.MS)},
                "group1", "group2",
                3,
                new List<double> { (756019.37596012) / ((23974.23487 + 807934.23874 + 23974.23487)/3), (234667.023876) / ((23974.23487 + 807934.23874 + 23974.23487)/3) },
                new List<double> { (38947.830453) / ((89754.37958 + 357698.893275 + 32587.4569384)/3), (8907235.0893645) / ((89754.37958 + 357698.893275 + 32587.4569384)/3) }
            }
        };
        [Test]
        [TestCaseSource("_calcStoichiometryBaselineCase")]
        public void PairwiseCompairison_calcStoichiometryBaselineCase_Pass(Peptide pep, List<Intensity> baseline, string g1, string g2, int minNumStoichiometries, List<double> stoich1, List<double> stoich2)
        {
            PairwiseCompairison PairwiseCompairisonTest = new PairwiseCompairison(pep, baseline, g1, g2, minNumStoichiometries);
            Assert.AreEqual(stoich1.Count(), PairwiseCompairisonTest.PeptideStoichiometriesGroupOne.Count());
            Assert.AreEqual(stoich2.Count(), PairwiseCompairisonTest.PeptideStoichiometriesGroupTwo.Count());
            for (int i=0; i < stoich1.Count(); i++)
            {
                Assert.That(stoich1[i], Is.EqualTo(PairwiseCompairisonTest.PeptideStoichiometriesGroupOne.Select(p => p.StoichiometryVal).ToList()[i]).Within(0.001));
            }
            for (int i = 0; i < stoich2.Count(); i++)
            {
                Assert.That(stoich2[i], Is.EqualTo(PairwiseCompairisonTest.PeptideStoichiometriesGroupTwo.Select(p => p.StoichiometryVal).ToList()[i]).Within(0.001));
            }
        }


        private static readonly object[] _calcStoichiometryPeptidePeptideCase =
        {
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file", "group1", 0, DetectionMS.MSMSIdentifiedButNotQuantified),
                        new Intensity("file", "group2", 500, DetectionMS.MS) },
                    new List<string>() { "group1", "group2" }, 1),
               new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file", "group1", 1000, DetectionMS.MSMS),
                        new Intensity("file", "group2", 500, DetectionMS.MS) },
                    new List<string>() { "group1", "group2" }, 1),
                "group1", "group2",
                3,
                new List<double> { },
                new List<double> { 1 }
            },
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file", "group1", 100, DetectionMS.MS),
                        new Intensity("file", "group2", 10, DetectionMS.MS) },
                    new List<string>() { "group1", "group2" }, 1),
               new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file", "group1", 1000, DetectionMS.MSMS),
                        new Intensity("file", "group2", 0, DetectionMS.NotDetected) },
                    new List<string>() { "group1", "group2" }, 1),
                "group1", "group2",
                3,
                new List<double> { 0.1 },
                new List<double> { }
            },
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file", "group1", 0, DetectionMS.MSMSIdentifiedButNotQuantified),
                        new Intensity("file", "group2", 500, DetectionMS.MS) },
                    new List<string>() { "group1", "group2" }, 1),
               new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file", "group1", 0, DetectionMS.NotDetected),
                        new Intensity("file", "group2", 500, DetectionMS.MS) },
                    new List<string>() { "group1", "group2" }, 1),
                "group1", "group2",
                3,
                new List<double> { },
                new List<double> { 1 }
            },
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file", "group1", 500, DetectionMS.MS),
                        new Intensity("file", "group2", 500, DetectionMS.MS) },
                    new List<string>() { "group1", "group2" }, 1),
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file", "group1", 1000, DetectionMS.MSMS),
                        new Intensity("file", "group2", 250, DetectionMS.MS) },
                    new List<string>() { "group1", "group2" }, 1),
                "group1", "group2",
                3,
                new List<double> { 0.5 },
                new List<double> { 2 }
            },
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file", "group1", 756019.37596012, DetectionMS.MS),
                        new Intensity("file", "group2", 234667.023876, DetectionMS.MS) },
                    new List<string>() { "group1", "group2" }, 1),
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file", "group1", 23974.23487, DetectionMS.MSMS),
                        new Intensity("file", "group2", 807934.23874, DetectionMS.MS) },
                    new List<string>() { "group1", "group2" }, 1),
                "group1", "group2",
                3,
                new List<double> { 756019.37596012/23974.23487 },
                new List<double> { 234667.023876/807934.23874 }
            },
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file1", "group1", 756019.37596012, DetectionMS.MS),
                        new Intensity("file2", "group1", 234667.023876, DetectionMS.MS),
                        new Intensity("file1", "group2", 38947.830453, DetectionMS.MS),
                        new Intensity("file2", "group2", 8907235.0893645, DetectionMS.MS),
                        new Intensity("file3", "group2", 0, DetectionMS.MSMSIdentifiedButNotQuantified),
                        new Intensity("file4", "group2", 8907235.0893645, DetectionMS.MS),
                        new Intensity("file1", "group3", 32974.039285, DetectionMS.MSMS),
                        new Intensity("file2", "group3", 8907235.0893645, DetectionMS.MS)
                    },
                    new List<string>() { "group1", "group2" }, 1),
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file1", "group1", 83475.98537, DetectionMS.MS),
                        new Intensity("file2", "group1", 3597.9823571, DetectionMS.MS),
                        new Intensity("file1", "group2", 97435.938475, DetectionMS.MS),
                        new Intensity("file2", "group2", 89745.43956, DetectionMS.MS),
                        new Intensity("file3", "group2", 398470.0974, DetectionMS.MSMS),
                        new Intensity("file4", "group2", 0, DetectionMS.NotDetected),
                        new Intensity("file1", "group3", 398470.0974, DetectionMS.MSMS),
                        new Intensity("file2", "group3", 23597.3659, DetectionMS.MS)
                    },
                    new List<string>() { "group1", "group2" }, 1),
                "group1", "group2", 3,
                new List<double> { 756019.37596012 / 83475.98537, 234667.023876 / 3597.9823571 },
                new List<double> { 38947.830453 / 97435.938475, 8907235.0893645 / 89745.43956 }
            }
        };
        [Test]
        [TestCaseSource("_calcStoichiometryPeptidePeptideCase")]
        public void PairwiseCompairison_calcStoichiometryPeptidePeptideCase_Pass(Peptide pep1, Peptide pep2, string g1, string g2, int minNumStoichiometries, List<double> stoich1, List<double> stoich2)
        {
            PairwiseCompairison PairwiseCompairisonTest = new PairwiseCompairison(pep1, pep2, g1, g2, minNumStoichiometries);
            Assert.AreEqual(stoich1.Count(), PairwiseCompairisonTest.PeptideStoichiometriesGroupOne.Count());
            Assert.AreEqual(stoich2.Count(), PairwiseCompairisonTest.PeptideStoichiometriesGroupTwo.Count());
            for (int i = 0; i < stoich1.Count(); i++)
            {
                Assert.That(stoich1[i], Is.EqualTo(PairwiseCompairisonTest.PeptideStoichiometriesGroupOne.Select(p => p.StoichiometryVal).ToList()[i]).Within(0.001));
            }
            for (int i = 0; i < stoich2.Count(); i++)
            {
                Assert.That(stoich2[i], Is.EqualTo(PairwiseCompairisonTest.PeptideStoichiometriesGroupTwo.Select(p => p.StoichiometryVal).ToList()[i]).Within(0.001));
            }
        }

        private static readonly object[] _calcMWStatsBaselineCase =
        {
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file", "group1", 50, DetectionMS.MS),
                        new Intensity("file", "group1", 50, DetectionMS.MS),
                        new Intensity("file", "group1", 50, DetectionMS.MS),
                        new Intensity("file", "group2", 100, DetectionMS.MS),
                        new Intensity("file", "group2", 100, DetectionMS.MS),
                        new Intensity("file", "group2", 100, DetectionMS.MS)
                    },
                    new List<string>() { "group1", "group2" }, 1),
                new List<Intensity> {
                new Intensity("file", "group1", 500, DetectionMS.MS),
                new Intensity("file", "group1", 500, DetectionMS.MS),
                new Intensity("file", "group1", 500, DetectionMS.MS),
                new Intensity("file", "group2", 10, DetectionMS.MS),
                new Intensity("file", "group2", 10, DetectionMS.MS),
                new Intensity("file", "group2", 10, DetectionMS.MS)},
                "group1", "group2", 3,
                0, 0.04685, 0.1, 10, 0.1, 10, 0.1, 10
            },
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 50, DetectionMS.MS),
                            new Intensity("file", "group1", 25, DetectionMS.MS),
                            new Intensity("file", "group1", 75, DetectionMS.MS),
                            new Intensity("file", "group1", 40, DetectionMS.MS),
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 75, DetectionMS.MS),
                            new Intensity("file", "group2", 100, DetectionMS.MS),
                            new Intensity("file", "group2", 400, DetectionMS.MS),
                            new Intensity("file", "group2", 800, DetectionMS.MS),
                            new Intensity("file", "group2", 600, DetectionMS.MS),
                            new Intensity("file", "group2", 300, DetectionMS.MS),
                            new Intensity("file", "group2", 200, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                    new List<Intensity> {
                    new Intensity("file", "group1", 500, DetectionMS.MS),
                    new Intensity("file", "group1", 450, DetectionMS.MS),
                    new Intensity("file", "group1", 250, DetectionMS.MS),
                    new Intensity("file", "group1", 300, DetectionMS.MS),
                    new Intensity("file", "group2", 10, DetectionMS.MS),
                    new Intensity("file", "group2", 20, DetectionMS.MS),
                    new Intensity("file", "group2", 50, DetectionMS.MS),
                    new Intensity("file", "group2", 40, DetectionMS.MS)},
                    "group1", "group2", 3,
                    0, 0.004998, 0.12, 11.66666, 0.053333, 3.33333333, 0.2, 26.66666
            },
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 3785016.398056, DetectionMS.MS),
                            new Intensity("file", "group1", 3783564.043746, DetectionMS.MS),
                            new Intensity("file", "group1", 984735.098, DetectionMS.MSMS),
                            new Intensity("file", "group1", 79854.09475, DetectionMS.MS),
                            new Intensity("file", "group1", 854875.9456, DetectionMS.MS),
                            new Intensity("file", "group1", 987354.098, DetectionMS.MSMS),
                            new Intensity("file", "group1", 398754.734, DetectionMS.MSMS),
                            new Intensity("file", "group1", 5768745.98735, DetectionMS.MS),
                            new Intensity("file", "group1", 987435.09846, DetectionMS.MS),
                            new Intensity("file", "group1", 847543.098754, DetectionMS.MSMS),

                            new Intensity("file", "group2", 9856.098436, DetectionMS.MS),
                            new Intensity("file", "group2", 6042.09476, DetectionMS.MS),
                            new Intensity("file", "group2", 96539.934085, DetectionMS.MS),
                            new Intensity("file", "group2", 9856.09865, DetectionMS.MS),
                            new Intensity("file", "group2", 9856.09756, DetectionMS.MS),
                            new Intensity("file", "group2", 48949.5476, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                    new List<Intensity> {
                    new Intensity("file", "group1", 897694.098743, DetectionMS.MS),
                    new Intensity("file", "group1", 24867.08197324, DetectionMS.MS),
                    new Intensity("file", "group1", 34875.0486, DetectionMS.MS),
                    new Intensity("file", "group1", 238475.340875, DetectionMS.MS),
                    new Intensity("file", "group1", 3434895.19087, DetectionMS.MS),
                    new Intensity("file", "group1", 45896.308475, DetectionMS.MS),
                    new Intensity("file", "group1", 408795.845, DetectionMS.MS),
                    new Intensity("file", "group1", 780243.087945, DetectionMS.MS),
                    new Intensity("file", "group1", 80793.98745, DetectionMS.MS),
                    new Intensity("file", "group1", 89743.987054, DetectionMS.MS),
                    new Intensity("file", "group1", 907384.09734, DetectionMS.MS),
                    new Intensity("file", "group1", 908743.34987, DetectionMS.MS),
                    new Intensity("file", "group1", 9876.98743, DetectionMS.MS),
                    new Intensity("file", "group1", 976543.09745, DetectionMS.MS),
                    new Intensity("file", "group1", 98073250.907435, DetectionMS.MS),
                    new Intensity("file", "group1", 97549.98754, DetectionMS.MS),
                    new Intensity("file", "group1", 987549.98754, DetectionMS.MS),
                    new Intensity("file", "group1", 865843.98463, DetectionMS.MS),

                    new Intensity("file", "group2", 896707.753743, DetectionMS.MS),
                    new Intensity("file", "group2", 23880.73697324, DetectionMS.MS),
                    new Intensity("file", "group2", 33888.7036, DetectionMS.MS),
                    new Intensity("file", "group2", 237488.995875, DetectionMS.MS),
                    new Intensity("file", "group2", 3433908.84587, DetectionMS.MS),
                    new Intensity("file", "group2", 44909.963475, DetectionMS.MS),
                    new Intensity("file", "group2", 407809.5, DetectionMS.MS),
                    new Intensity("file", "group2", 779256.742945, DetectionMS.MS),
                    new Intensity("file", "group2", 79807.64245, DetectionMS.MS),
                    new Intensity("file", "group2", 88757.642054, DetectionMS.MS),
                    new Intensity("file", "group2", 906397.75234, DetectionMS.MS),
                    new Intensity("file", "group2", 907757.00487, DetectionMS.MS),
                    new Intensity("file", "group2", 8890.64243, DetectionMS.MS),
                    new Intensity("file", "group2", 975556.75245, DetectionMS.MS),
                    new Intensity("file", "group2", 98072264.562435, DetectionMS.MS),
                    new Intensity("file", "group2", 96563.64254, DetectionMS.MS),
                    new Intensity("file", "group2", 986563.64254, DetectionMS.MS),
                    new Intensity("file", "group2", 864857.63963, DetectionMS.MS)},
                    "group1", "group2", 3,
                    59, 0.0004995, 0.163037939, 0.001629926, 0.013203507, 0.000999196, 0.953835614, 0.015965038
            }

        };
        [Test]
        [TestCaseSource("_calcMWStatsBaselineCase")]
        public void PairwiseCompairison_calcMWStatsBaselineCase_Pass(Peptide pep, List<Intensity> baseline, string g1, string g2, int minNumStoichiometries, 
            double mwStat, double mwPVal, double g1Med, double g2Med, double g1Min, double g2Min, double g1Max, double g2Max)
        {
            PairwiseCompairison PairwiseCompairisonTest = new PairwiseCompairison(pep, baseline, g1, g2, minNumStoichiometries);
            Assert.That(mwStat, Is.EqualTo(PairwiseCompairisonTest.MWStat).Within(0.01));
            Assert.That(mwPVal, Is.EqualTo(PairwiseCompairisonTest.MWPVal).Within(0.01));
            Assert.That(g1Med, Is.EqualTo(PairwiseCompairisonTest.PeptideStoichiometriesGroupOneMedian).Within(0.001));
            Assert.That(g2Med, Is.EqualTo(PairwiseCompairisonTest.PeptideStoichiometriesGroupTwoMedian).Within(0.001));
            Assert.That(g1Min, Is.EqualTo(PairwiseCompairisonTest.PeptideStoichiometriesGroupOneMin).Within(0.001));
            Assert.That(g2Min, Is.EqualTo(PairwiseCompairisonTest.PeptideStoichiometriesGroupTwoMin).Within(0.001));
            Assert.That(g1Max, Is.EqualTo(PairwiseCompairisonTest.PeptideStoichiometriesGroupOneMax).Within(0.001));
            Assert.That(g2Max, Is.EqualTo(PairwiseCompairisonTest.PeptideStoichiometriesGroupTwoMax).Within(0.001));
        }

    }
}
