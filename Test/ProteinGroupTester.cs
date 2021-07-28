using NUnit.Framework;
using System;
using System.Collections.Generic;
using PTMStoichiometry;
using System.Linq;
using System.Text;
using System.IO;

namespace Test
{
    [TestFixture]
    class ProteinGroupTester
    {

        

        private static readonly object[] _getBaseLinePeptides =
        {
            new object[] {
                "Prot",
                new List<Peptide>
                {
                    new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 10),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 15),
                            new Intensity("file", "group1", 30),

                            new Intensity("file", "group2", 10),
                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 15),
                            new Intensity("file", "group2", 30) 
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq2", "Seq2", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 30),
                            new Intensity("file", "group1", 60),

                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 30),
                            new Intensity("file", "group2", 60) 
                        },
                    new List<string>() { "group1", "group2" }, 1)
                },
                new List<string> { "group1", "group2" },
                2, 0, 2, true, 2, 0.5, false, 3, "group1",
                new List<string> { "Seq1", "Seq2" }
            },
            new object[] {
                "Prot",
                new List<Peptide>
                {
                    new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 10),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 15),
                            new Intensity("file", "group1", 30),

                            new Intensity("file", "group2", 10),
                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 15),
                            new Intensity("file", "group2", 30)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq2", "Seq2", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 30),
                            new Intensity("file", "group1", 60),

                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 30),
                            new Intensity("file", "group2", 60)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq3", "Seq3", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 100),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 50),
                            new Intensity("file", "group1", 600),

                            new Intensity("file", "group2", 80),
                            new Intensity("file", "group2", 600),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 3)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq4", "Seq4", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 30),
                            new Intensity("file", "group1", 60),

                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 30),
                            new Intensity("file", "group2", 60)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq4", "Seq4", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 30),
                            new Intensity("file", "group1", 60),

                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 30),
                            new Intensity("file", "group2", 60)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq5a", "Seq5", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 30),
                            new Intensity("file", "group1", 60),

                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 30),
                            new Intensity("file", "group2", 60)
                        },
                        new List<string>() { "group1", "group2" }, 1)
                },
                new List<string> { "group1", "group2" },
                2, 0, 2, true, 2, 0.5, false, 3, "group1",
                new List<string> { "Seq1", "Seq2" }
            },
            new object[] {
                "Prot",
                new List<Peptide>
                {
                    new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 10),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 15),
                            new Intensity("file", "group1", 30),

                            new Intensity("file", "group2", 10),
                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 15),
                            new Intensity("file", "group2", 30),

                            new Intensity("file", "group3", 10),
                            new Intensity("file", "group3", 20),
                            new Intensity("file", "group3", 15),
                            new Intensity("file", "group3", 30)
                        },
                        new List<string>() { "group1", "group2", "group3" }, 1),
                   new Peptide("Seq2", "Seq2", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 30),
                            new Intensity("file", "group1", 60),

                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 30),
                            new Intensity("file", "group2", 60),

                            new Intensity("file", "group3", 20),
                            new Intensity("file", "group3", 40),
                            new Intensity("file", "group3", 30),
                            new Intensity("file", "group3", 60)
                        },
                        new List<string>() { "group1", "group2" , "group3"}, 1),
                   new Peptide("Seq3", "Seq3", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 100),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 50),
                            new Intensity("file", "group1", 600),

                            new Intensity("file", "group2", 80),
                            new Intensity("file", "group2", 600),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 3),

                            new Intensity("file", "group3", 100),
                            new Intensity("file", "group3", 20),
                            new Intensity("file", "group3", 50),
                            new Intensity("file", "group3", 600)

                        }, 
                        new List<string>() { "group1", "group2" , "group3"}, 1)
                },
                new List<string> { "group1", "group2" , "group3"},
                2, 0, 2, true, 2, 0.5, false, 3, "group1",
                new List<string> { "Seq1", "Seq2" }
            },
            new object[] {
                "Prot",
                new List<Peptide>
                {
                    new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 10),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 15),
                            new Intensity("file", "group1", 30),

                            new Intensity("file", "group2", 10),
                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 15),
                            new Intensity("file", "group2", 30)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq2", "Seq2", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 30),
                            new Intensity("file", "group1", 60),

                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 30),
                            new Intensity("file", "group2", 60)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq3", "Seq3", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 100),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 50),
                            new Intensity("file", "group1", 600),

                            new Intensity("file", "group2", 80),
                            new Intensity("file", "group2", 600),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 3)
                        },
                        new List<string>() { "group1", "group2"}, 1),
                   new Peptide("Seq4", "Seq4", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 150),
                            new Intensity("file", "group1", 50),
                            new Intensity("file", "group1", 75),
                            new Intensity("file", "group1", 200),

                            new Intensity("file", "group2", 90),
                            new Intensity("file", "group2", 100),
                            new Intensity("file", "group2", 50),
                            new Intensity("file", "group2", 10)
                        },
                        new List<string>() { "group1", "group2"}, 1),
                   new Peptide("Seq5", "Seq5", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 160),
                            new Intensity("file", "group1", 55),
                            new Intensity("file", "group1", 80),
                            new Intensity("file", "group1", 240),

                            new Intensity("file", "group2", 110),
                            new Intensity("file", "group2", 120),
                            new Intensity("file", "group2", 75),
                            new Intensity("file", "group2", 24)
                        },
                        new List<string>() { "group1", "group2"}, 1)

                },
                new List<string> { "group1", "group2" },
                2, 0, 2, true, 2, 0.5, false, 3, "group1",
                new List<string> { "Seq3", "Seq4", "Seq5" }
            },
            new object[] {
                "Prot",
                new List<Peptide>
                {
                    new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 10),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 15),
                            new Intensity("file", "group1", 30),

                            new Intensity("file", "group2", 10),
                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 15),
                            new Intensity("file", "group2", 30)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq2", "Seq2", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 30),
                            new Intensity("file", "group1", 60),

                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 30),
                            new Intensity("file", "group2", 60)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq3", "Seq3", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 100),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 50),
                            new Intensity("file", "group1", 600),

                            new Intensity("file", "group2", 80),
                            new Intensity("file", "group2", 600),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 3)
                        },
                        new List<string>() { "group1", "group2"}, 1),
                   new Peptide("Seq4", "Seq4", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 150),
                            new Intensity("file", "group1", 50),
                            new Intensity("file", "group1", 75),
                            new Intensity("file", "group1", 200),

                            new Intensity("file", "group2", 90),
                            new Intensity("file", "group2", 100),
                            new Intensity("file", "group2", 50),
                            new Intensity("file", "group2", 10)
                        },
                        new List<string>() { "group1", "group2"}, 1),
                   new Peptide("Seq5", "Seq5", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 160),
                            new Intensity("file", "group1", 55),
                            new Intensity("file", "group1", 80),
                            new Intensity("file", "group1", 240),

                            new Intensity("file", "group2", 110),
                            new Intensity("file", "group2", 120),
                            new Intensity("file", "group2", 75),
                            new Intensity("file", "group2", 24)
                        },
                        new List<string>() { "group1", "group2"}, 1)

                },
                new List<string> { "group1", "group2" },
                2, 0, 2, true, 2, 0.7, false, 3, "group1",
                new List<string> { "Seq1", "Seq2" }
            }
        };

        [Test]
        [TestCaseSource("_getBaseLinePeptides")]
        public void ProteinGroup_getBaseLinePeptides_Pass(string proteinAccession, List<Peptide> peptides, List<string> groups, int reqNumUnmodPeptides, int reqNumModPeptides, int reqNumOfPepeptides,
            int reqNumBaselinePeptides, double correlationCutOff, bool compareUnmod, int minNumStoichiometries, string groupToCompare, List<string> baselinePeptideSeq)
        {
           // Extensions.IncludeSharedPeptides(peptides, true); //set isUnique
            ProteinGroup ProteinGroupAllGroupsTest = new ProteinGroup(proteinAccession, peptides, groups, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides, 
                reqNumBaselinePeptides, 3, correlationCutOff, compareUnmod, minNumStoichiometries);
            ProteinGroup ProteinGroupSetGroupTest = new ProteinGroup(proteinAccession, peptides, groups, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides,
                reqNumBaselinePeptides, 3, correlationCutOff, compareUnmod, minNumStoichiometries, groupToCompare);

            Assert.AreEqual(baselinePeptideSeq.Count(), ProteinGroupAllGroupsTest.BaselinePeptides.Count());
            Assert.AreEqual(baselinePeptideSeq.Count(), ProteinGroupSetGroupTest.BaselinePeptides.Count());
            Assert.AreEqual(baselinePeptideSeq, ProteinGroupAllGroupsTest.BaselinePeptides.Select(p => p.Sequence));
            Assert.AreEqual(baselinePeptideSeq, ProteinGroupSetGroupTest.BaselinePeptides.Select(p => p.Sequence));
        }

        private static readonly object[] _useProt =
        {
            new object[] {
                "Prot",
                new List<Peptide>
                {
                    new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 10),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 15),
                            new Intensity("file", "group1", 30),

                            new Intensity("file", "group2", 10),
                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 15),
                            new Intensity("file", "group2", 30)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq2", "Seq2", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 30),
                            new Intensity("file", "group1", 60),

                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 30),
                            new Intensity("file", "group2", 60)
                        },
                    new List<string>() { "group1", "group2" }, 1)
                },
                new List<string> { "group1", "group2" },
                3, 1, 4, 0.5, false, 3, "group1",
                false
            },
            new object[] {
                "Prot",
                new List<Peptide>
                {
                    new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 10),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 15),
                            new Intensity("file", "group1", 30),

                            new Intensity("file", "group2", 10),
                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 15),
                            new Intensity("file", "group2", 30)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq2", "Seq2", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 30),
                            new Intensity("file", "group1", 60),

                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 30),
                            new Intensity("file", "group2", 60)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq3", "Seq3", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 100),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 50),
                            new Intensity("file", "group1", 600),

                            new Intensity("file", "group2", 80),
                            new Intensity("file", "group2", 600),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 3)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq4", "Seq4", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 30),
                            new Intensity("file", "group1", 60),

                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 30),
                            new Intensity("file", "group2", 60)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq4", "Seq4", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 30),
                            new Intensity("file", "group1", 60),

                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 30),
                            new Intensity("file", "group2", 60)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq5a", "Seq5", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 30),
                            new Intensity("file", "group1", 60),

                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 30),
                            new Intensity("file", "group2", 60)
                        },
                        new List<string>() { "group1", "group2" }, 1)
                },
                new List<string> { "group1", "group2" },
                3, 1, 4, 0.5, false, 3, "group1",
                false
            },
            new object[] {
                "Prot",
                new List<Peptide>
                {
                    new Peptide("Seq1", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 10),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 15),
                            new Intensity("file", "group1", 30),

                            new Intensity("file", "group2", 10),
                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 15),
                            new Intensity("file", "group2", 30)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq2", "Seq2", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 30),
                            new Intensity("file", "group1", 60),

                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 30),
                            new Intensity("file", "group2", 60)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq3", "Seq3", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 100),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 50),
                            new Intensity("file", "group1", 600),

                            new Intensity("file", "group2", 80),
                            new Intensity("file", "group2", 600),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 3)
                        },
                        new List<string>() { "group1", "group2"}, 1),
                   new Peptide("Seq4", "Seq4", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 150),
                            new Intensity("file", "group1", 50),
                            new Intensity("file", "group1", 75),
                            new Intensity("file", "group1", 200),

                            new Intensity("file", "group2", 90),
                            new Intensity("file", "group2", 100),
                            new Intensity("file", "group2", 50),
                            new Intensity("file", "group2", 10)
                        },
                        new List<string>() { "group1", "group2"}, 1),
                   new Peptide("Seq5a", "Seq5", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 160),
                            new Intensity("file", "group1", 55),
                            new Intensity("file", "group1", 80),
                            new Intensity("file", "group1", 240),

                            new Intensity("file", "group2", 110),
                            new Intensity("file", "group2", 120),
                            new Intensity("file", "group2", 75),
                            new Intensity("file", "group2", 24)
                        },
                        new List<string>() { "group1", "group2"}, 1)

                },
                new List<string> { "group1", "group2" },
                3, 1, 4, 0.5, false, 3, "group1",
                false
            }
        };



        [Test]
        [TestCaseSource("_useProt")]
        public void ProteinGroup_useProt_Pass(string proteinAccession, List<Peptide> peptides, List<string> groups, int reqNumUnmodPeptides, int reqNumModPeptides, int reqNumOfPepeptides,
            double correlationCutOff, bool compareUnmod, int minNumStoichiometries, string groupToCompare, bool useProtBaseline)
        {
            ProteinGroup ProteinGroupBaselineAllGroupTest = new ProteinGroup(proteinAccession, peptides, groups, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides,
                reqNumUnmodPeptides, 3, correlationCutOff, compareUnmod, minNumStoichiometries);
            ProteinGroup ProteinGroupBaselineSetGroupTest = new ProteinGroup(proteinAccession, peptides, groups, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides,
                reqNumUnmodPeptides, 3, correlationCutOff, compareUnmod, minNumStoichiometries, groupToCompare);

            Assert.AreEqual(useProtBaseline, ProteinGroupBaselineAllGroupTest.useProt);
            Assert.AreEqual(useProtBaseline, ProteinGroupBaselineSetGroupTest.useProt);
        }



        private static readonly object[] _calCompairison =
 {
            new object[] {
                "Prot",
                new List<Peptide>
                {
                    new Peptide("Seq1a", "Seq1", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 10),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 15),
                            new Intensity("file", "group1", 30),

                            new Intensity("file", "group2", 10),
                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 15),
                            new Intensity("file", "group2", 30),

                            new Intensity("file", "group3", 10),
                            new Intensity("file", "group3", 20),
                            new Intensity("file", "group3", 15),
                            new Intensity("file", "group3", 30)
                        },
                        new List<string>() { "group1", "group2", "group3" }, 1),
                   new Peptide("Seq2b", "Seq2", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 30),
                            new Intensity("file", "group1", 60),

                            new Intensity("file", "group2", 20),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 30),
                            new Intensity("file", "group2", 60),

                            new Intensity("file", "group3", 20),
                            new Intensity("file", "group3", 40),
                            new Intensity("file", "group3", 30),
                            new Intensity("file", "group3", 60)
                        },
                        new List<string>() { "group1", "group2", "group3" }, 1),
                   new Peptide("Seq3", "Seq3", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 100),
                            new Intensity("file", "group1", 20),
                            new Intensity("file", "group1", 50),
                            new Intensity("file", "group1", 600),

                            new Intensity("file", "group2", 80),
                            new Intensity("file", "group2", 600),
                            new Intensity("file", "group2", 40),
                            new Intensity("file", "group2", 3),

                            new Intensity("file", "group3", 80),
                            new Intensity("file", "group3", 600),
                            new Intensity("file", "group3", 40),
                            new Intensity("file", "group3", 3)
                        },
                        new List<string>() { "group1", "group2", "group3" }, 1),
                   new Peptide("Seq4", "Seq4", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 150),
                            new Intensity("file", "group1", 50),
                            new Intensity("file", "group1", 75),
                            new Intensity("file", "group1", 200),

                            new Intensity("file", "group2", 90),
                            new Intensity("file", "group2", 100),
                            new Intensity("file", "group2", 50),
                            new Intensity("file", "group2", 10),

                            new Intensity("file", "group3", 90),
                            new Intensity("file", "group3", 100),
                            new Intensity("file", "group3", 50),
                            new Intensity("file", "group3", 10)
                        },
                        new List<string>() { "group1", "group2", "group3" }, 1),
                   new Peptide("Seq5", "Seq5", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 160),
                            new Intensity("file", "group1", 55),
                            new Intensity("file", "group1", 80),
                            new Intensity("file", "group1", 240),

                            new Intensity("file", "group2", 110),
                            new Intensity("file", "group2", 120),
                            new Intensity("file", "group2", 75),
                            new Intensity("file", "group2", 24),

                            new Intensity("file", "group3", 110),
                            new Intensity("file", "group3", 120),
                            new Intensity("file", "group3", 75),
                            new Intensity("file", "group3", 24)
                        },
                        new List<string>() { "group1", "group2", "group3" }, 1)

                },
                new List<string> { "group1", "group2", "group3" }, "group1",
                6, 4, 30, 20
            }
        };

        //calcComparison - since PairwiseCompairison is tested elsewhere, just want to check is returning correct number of PairwiseCompairisons depending on inputs 
        [Test]
        [TestCaseSource("_calCompairison")]
        public void ProteinGroup_calcComparison_Pass(string proteinAccession, List<Peptide> peptides, List<string> groups, string groupToCompare,
            int numPairwiseCompairisonsBaselineAllGroup, int numPairwiseCompairisonsBaselineSetGroup, int numPairwiseCompairisonsPeptidePeptideAllGroup,
            int numPairwiseCompairisonsPeptidePeptideSetGroup)
        {

            ProteinGroup ProteinGroupBaselineAllGroupTest = new ProteinGroup(proteinAccession, peptides, groups, 3, 1, 4,
                3, 3, 0.5, false, 3);
            ProteinGroup ProteinGroupBaselineSetGroupTest = new ProteinGroup(proteinAccession, peptides, groups, 3, 1, 4,
                3, 3, 0.5, false, 3, groupToCompare);

            Assert.AreEqual(numPairwiseCompairisonsBaselineAllGroup, ProteinGroupBaselineAllGroupTest.ProteinPairwiseComparisons.Count());
            Assert.AreEqual(numPairwiseCompairisonsBaselineSetGroup, ProteinGroupBaselineSetGroupTest.ProteinPairwiseComparisons.Count());

        }
    }
}

