using NUnit.Framework;
using System;
using System.Collections.Generic;
using PTMStoichiometry;
using System.Linq;

namespace Test
{
    [TestFixture]
    class PairwiseCompairisonTester
    {
        
        private static readonly object[] _calcStoichiometryBaselineCase =
        {
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file", "group1", 0),
                        new Intensity("file", "group2", 500) },
                    new List<string>() { "group1", "group2" }, 1),
                new List<Intensity> {
                new Intensity("file", "group1", 1000),
                new Intensity("file", "group2", 500)},
                "group1", "group2",
                3,
                new List<double> { 0 },
                new List<double> { 1 }
            },
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism", 
                    new List<Intensity>() { 
                        new Intensity("file", "group1", 500), 
                        new Intensity("file", "group2", 500) }, 
                    new List<string>() { "group1", "group2" }, 1),
                new List<Intensity> {
                new Intensity("file", "group1", 1000),
                new Intensity("file", "group2", 250)},
                "group1", "group2",
                3,
                new List<double> { 0.5 },
                new List<double> { 2 }
            },
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism", 
                    new List<Intensity>() { 
                        new Intensity("file", "group1", 756019.37596012), 
                        new Intensity("file", "group2", 234667.023876) },
                    new List<string>() { "group1", "group2" }, 1),
                new List<Intensity> {
                new Intensity("file", "group1", 23974.23487),
                new Intensity("file", "group2", 807934.23874)},
                "group1", "group2",
                3,
                new List<double> { 756019.37596012/23974.23487 },
                new List<double> { 234667.023876/807934.23874 }
            },
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism", 
                    new List<Intensity>() { 
                        new Intensity("file", "group1", 756019.37596012), 
                        new Intensity("file", "group1", 234667.023876),
                        new Intensity("file", "group2", 38947.830453),
                        new Intensity("file", "group2", 8907235.0893645),
                        new Intensity("file", "group2", 0),
                        new Intensity("file", "group3", 8907235.0893645),
                        new Intensity("file", "group3", 3409750.2309)
                    },
                    new List<string>() { "group1", "group2" }, 1),
                new List<Intensity> {
                new Intensity("file", "group1", 23974.23487),
                new Intensity("file", "group1", 807934.23874),
                new Intensity("file", "group1", 23974.23487),
                new Intensity("file", "group2", 89754.37958),
                new Intensity("file", "group2", 357698.893275),
                new Intensity("file", "group2", 32587.4569384),
                new Intensity("file", "group3", 89754.37958),
                new Intensity("file", "group3", 357698.893275),
                new Intensity("file", "group3", 32587.4569384)},
                "group1", "group2",
                3,
                new List<double> { (756019.37596012) / 23974.23487, (234667.023876) / 23974.23487 },
                new List<double> { (38947.830453) / 89754.37958, (8907235.0893645) / 89754.37958, 0 }
            }
        };
        /// <summary>
        /// Test that checks that calcStoichiometry in the baseline case is returning the correct stoichiometry values
        /// </summary>
        /// <param name="pep">peptide to compare (numerator)</param>
        /// <param name="baseline">the baseline value (denominator)</param>
        /// <param name="g1">group to be compared</param>
        /// <param name="g2">group to be compared</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        /// <param name="stoich1">correct stoichiometries in g1</param>
        /// <param name="stoich2">correct stoichiometries in g2</param>
        [Test]
        [TestCaseSource("_calcStoichiometryBaselineCase")]
        public void PairwiseCompairison_calcStoichiometryBaselineCase_Pass(Peptide pep, List<Intensity> baseline, string g1, string g2, 
            int minNumStoichiometries, List<double> stoich1, List<double> stoich2)
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

  
        private static readonly object[] _calcStoichiometryPTMBaselineCase =
        {
            new object[] {
                new List<Peptide> 
                {
                    new Peptide("Seq1moda", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file", "group1", 100),
                        new Intensity("file1", "group2", 500) },
                    new List<string>() { "group1", "group2" }, 1)
                },
                new List<Intensity> {
                new Intensity("file", "group1", 1000),
                new Intensity("file1", "group2", 500)},
                "group1", "group2", "moda",
                3,
                new List<double> { 0.1 },
                new List<double> { 1 }
                
            },

            new object[] {
                new List<Peptide>
                {
                    new Peptide("Seq1moda", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 500),
                            new Intensity("file2", "group2", 500) },
                        new List<string>() { "group1", "group2" }, 1)
                },
                new List<Intensity> {
                new Intensity("file", "group1", 1000),
                new Intensity("file2", "group2", 250)},
                "group1", "group2", "moda",
                3,
                new List<double> { 0.5 },
                new List<double> { 2 }
            },


            new object[] {
                new List<Peptide>
                {
                    new Peptide("Seq1modb", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file1", "group1", 500),
                            new Intensity("file2", "group2", 200),
                            new Intensity("file11", "group1", 500),
                            new Intensity("file22", "group2", 200)},
                        new List<string>() { "group1", "group2" }, 1),
                    new Peptide("Seq1modb", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file1", "group1", 400),
                            new Intensity("file22", "group2", 400),
                            new Intensity("file11", "group1", 0),
                            new Intensity("file2", "group2", 0)},
                        new List<string>() { "group1", "group2" }, 1),
                    new Peptide("Seq1modb", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file11", "group1", 100),
                            new Intensity("file2", "group2", 600),
                            new Intensity("file1", "group1", 0),
                            new Intensity("file22", "group2", 0)},
                        new List<string>() { "group1", "group2" }, 1)
                },
                new List<Intensity> {
                new Intensity("file1", "group1", 1000),
                new Intensity("file11", "group1", 1000),
                new Intensity("file2", "group2", 400),
                new Intensity("file22", "group2", 400)},
                "group1", "group2", "modb",
                3,
                new List<double> { 0.9, 0.6 },
                new List<double> { 2, 1.5 }
            }
        };
        /// <summary>
        /// Test that checks that calcStoichiometry in the baseline case is returning the correct stoichiometry values for PTMs
        /// </summary>
        /// <param name="peps">peptides to compare (numerator)</param>
        /// <param name="baseline">the baseline value (denominator)</param>
        /// <param name="g1">group to be compared</param>
        /// <param name="g2">group to be compared</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        /// <param name="stoich1">correct stoichiometries in g1</param>
        /// <param name="stoich2">correct stoichiometries in g2</param>
        [Test]
        [TestCaseSource("_calcStoichiometryPTMBaselineCase")]
        public void PairwiseCompairison_calcStoichiometryPTMBaselineCase_Pass(List<Peptide> peps, List<Intensity> baseline, string g1, string g2,
            string ptm, int minNumStoichiometries, List<double> stoich1, List<double> stoich2)
        {
            
            PairwiseCompairison PairwiseCompairisonTest = new PairwiseCompairison(peps, baseline, g1, g2, minNumStoichiometries, ptm);
            Assert.AreEqual(stoich1.Count(), PairwiseCompairisonTest.PeptideStoichiometriesGroupOne.Count());
            Assert.AreEqual(stoich2.Count(), PairwiseCompairisonTest.PeptideStoichiometriesGroupTwo.Count());
            for (int i = 0; i < stoich1.Count(); i++)
            {
                Assert.That(PairwiseCompairisonTest.PeptideStoichiometriesGroupOne.Select(p => p.StoichiometryVal).ToList()[i], Is.EqualTo(stoich1[i]).Within(0.001));
            }
            for (int i = 0; i < stoich2.Count(); i++)
            {
                Assert.That(PairwiseCompairisonTest.PeptideStoichiometriesGroupTwo.Select(p => p.StoichiometryVal).ToList()[i], Is.EqualTo(stoich2[i]).Within(0.001));
            }
        }


       

        private static readonly object[] _calcMWStatsBaselineCase =
        {
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                    new List<Intensity>() {
                        new Intensity("file", "group1", 50),
                        new Intensity("file", "group1", 50),
                        new Intensity("file", "group1", 50),
                        new Intensity("file", "group2", 100),
                        new Intensity("file", "group2", 100),
                        new Intensity("file", "group2", 100)
                    },
                    new List<string>() { "group1", "group2" }, 1),
                new List<Intensity> {
                new Intensity("file", "group1", 500),
                new Intensity("file", "group1", 500),
                new Intensity("file", "group1", 500),
                new Intensity("file", "group2", 10),
                new Intensity("file", "group2", 10),
                new Intensity("file", "group2", 10)},
                "group1", "group2", 3,
                0, 0.04685, 0.1, 10, 0.1, 10, 0.1, 10
            },
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 50),
                            new Intensity("file", "group1", 25),
                            new Intensity("file", "group1", 75),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 75),
                            new Intensity("file", "group2", 100),
                            new Intensity("file", "group2", 400),
                            new Intensity("file", "group2", 800),
                            new Intensity("file", "group2", 600),
                            new Intensity("file", "group2", 300),
                            new Intensity("file", "group2", 200)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                    new List<Intensity> {
                    new Intensity("file", "group1", 500),
                    new Intensity("file", "group1", 450),
                    new Intensity("file", "group1", 250),
                    new Intensity("file", "group1", 300),
                    new Intensity("file", "group2", 10),
                    new Intensity("file", "group2", 20),
                    new Intensity("file", "group2", 50),
                    new Intensity("file", "group2", 40)},
                    "group1", "group2", 3,
                    0, 0.004998, 0.12, 11.66666, 0.053333, 3.33333333, 0.2, 26.66666
            },
            new object[] {
                new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 3785016.398056),
                            new Intensity("file", "group1", 3783564.043746),
                            new Intensity("file", "group1", 984735.098),
                            new Intensity("file", "group1", 79854.09475),
                            new Intensity("file", "group1", 854875.9456),
                            new Intensity("file", "group1", 987354.098),
                            new Intensity("file", "group1", 398754.734),
                            new Intensity("file", "group1", 5768745.98735),
                            new Intensity("file", "group1", 987435.09846),
                            new Intensity("file", "group1", 847543.098754),

                            new Intensity("file", "group2", 9856.098436),
                            new Intensity("file", "group2", 6042.09476),
                            new Intensity("file", "group2", 96539.934085),
                            new Intensity("file", "group2", 9856.09865),
                            new Intensity("file", "group2", 9856.09756),
                            new Intensity("file", "group2", 48949.5476)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                    new List<Intensity> {
                    new Intensity("file", "group1", 897694.098743),
                    new Intensity("file", "group1", 24867.08197324),
                    new Intensity("file", "group1", 34875.0486),
                    new Intensity("file", "group1", 238475.340875),
                    new Intensity("file", "group1", 3434895.19087),
                    new Intensity("file", "group1", 45896.308475),
                    new Intensity("file", "group1", 408795.845),
                    new Intensity("file", "group1", 780243.087945),
                    new Intensity("file", "group1", 80793.98745),
                    new Intensity("file", "group1", 89743.987054),
                    new Intensity("file", "group1", 907384.09734),
                    new Intensity("file", "group1", 908743.34987),
                    new Intensity("file", "group1", 9876.98743),
                    new Intensity("file", "group1", 976543.09745),
                    new Intensity("file", "group1", 98073250.907435),
                    new Intensity("file", "group1", 97549.98754),
                    new Intensity("file", "group1", 987549.98754),
                    new Intensity("file", "group1", 865843.98463),

                    new Intensity("file", "group2", 896707.753743),
                    new Intensity("file", "group2", 23880.73697324),
                    new Intensity("file", "group2", 33888.7036),
                    new Intensity("file", "group2", 237488.995875),
                    new Intensity("file", "group2", 3433908.84587),
                    new Intensity("file", "group2", 44909.963475),
                    new Intensity("file", "group2", 407809.5),
                    new Intensity("file", "group2", 779256.742945),
                    new Intensity("file", "group2", 79807.64245),
                    new Intensity("file", "group2", 88757.642054),
                    new Intensity("file", "group2", 906397.75234),
                    new Intensity("file", "group2", 907757.00487),
                    new Intensity("file", "group2", 8890.64243),
                    new Intensity("file", "group2", 975556.75245),
                    new Intensity("file", "group2", 98072264.562435),
                    new Intensity("file", "group2", 96563.64254),
                    new Intensity("file", "group2", 986563.64254),
                    new Intensity("file", "group2", 864857.63963)},
                    "group1", "group2", 3,
                    59, 0.0004995, 1.658557295, 0.016605811, 0.13431704, 0.010179878, 9.703207906, 0.162652985
            }

        };
        /// <summary>
        /// Test to check that calcMWStats is returning the correct values
        /// </summary>
        /// <param name="pep">peptide to compare (numerator)</param>
        /// <param name="baseline">baseline (denominator)</param>
        /// <param name="g1">group to be compared</param>
        /// <param name="g2">group to be compared</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        /// <param name="mwStat">Mann-Whitney statistic (calculated in R)</param>
        /// <param name="mwPVal">Mann-Whitney p-value (calculated in R)</param>
        /// <param name="g1Med">median of g1 stoichiometries</param>
        /// <param name="g2Med">median of g2 stoichiometries</param>
        /// <param name="g1Min">minimum of g1 stoichiometries</param>
        /// <param name="g2Min">minimum of g2 stoichiometries</param>
        /// <param name="g1Max">maximum of g1 stoichiometries</param>
        /// <param name="g2Max">maximum of g2 stoichiometries</param>
        [Test]
        [TestCaseSource("_calcMWStatsBaselineCase")]
        public void PairwiseCompairison_calcMWStatsBaselineCase_Pass(Peptide pep, List<Intensity> baseline, string g1, string g2, int minNumStoichiometries, 
            double mwStat, double mwPVal, double g1Med, double g2Med, double g1Min, double g2Min, double g1Max, double g2Max)
        {
            PairwiseCompairison PairwiseCompairisonTest = new PairwiseCompairison(pep, baseline, g1, g2, minNumStoichiometries);
            Assert.That(mwStat, Is.EqualTo(PairwiseCompairisonTest.MWStat).Within(0.001));
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
