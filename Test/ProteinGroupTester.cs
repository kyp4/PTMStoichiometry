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
                2, 0, 2, 2, 0.5, false, 3, "group1",
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
                            new Intensity("file", "group1", 0.2),
                            new Intensity("file", "group1", 40),
                            new Intensity("file", "group1", 300),
                            new Intensity("file", "group1", 70),

                            new Intensity("file", "group2", 70),
                            new Intensity("file", "group2", 70),
                            new Intensity("file", "group2", 70),
                            new Intensity("file", "group2", 70)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq5[Acetylation: on X]", "Seq5", "Prot", "Gene", "Organism",
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
                1, 1, 2, 2, 0.5, false, 3, "group1",
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
                2, 0, 2, 2, 0.5, false, 3, "group1",
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
                2, 0, 2, 2, 0.5, false, 3, "group1",
                new List<string> { "Seq3", "Seq4", "Seq5" }
            }
            
        };
        /// <summary>
        /// Test to check that the correct baseline peptides are being selected
        /// </summary>
        /// <param name="proteinAccession">the protein acession</param>
        /// <param name="peptides">list of all peptides</param>
        /// <param name="groups">list of all groups</param>
        /// <param name="reqNumUnmodPeptides">min num of unmodified peptides that must be observed for a protein in order to consider the protein group (default=1)</param>
        /// <param name="reqNumModPeptides">min num of modified peptides that must be observed for a protein in order to consider the protein group (default=3)</param>
        /// <param name="reqNumOfPepeptides">min num of peptides that must be observed for a protein in order to consider the protein group </param>
        /// <param name="reqNumBaselinePeptides">min num of baseline peptides that must be observed for a protein in order to consider it (default=3)</param>
        /// <param name="correlationCutOff">min value at which two peptides will be considered to be correlated</param>
        /// <param name="compareUnmod">if false (default) only compare modified peptides to baseline, not unmodified peptides</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        /// <param name="groupToCompare">single group to compare all other groups to (not a required parameter)</param>
        /// <param name="baselinePeptideSeq">list of correct baseline peptide sequences to compare against</param>
        [Test]
        [TestCaseSource("_getBaseLinePeptides")]
        public void ProteinGroup_getBaseLinePeptides_Pass(string proteinAccession, List<Peptide> peptides, List<string> groups, int reqNumUnmodPeptides,
            int reqNumModPeptides, int reqNumOfPepeptides, int reqNumBaselinePeptides, double correlationCutOff, bool compareUnmod, int minNumStoichiometries,
            string groupToCompare, List<string> baselinePeptideSeq)
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


        /// <summary>
        /// Test to check that useProt is being set correctly
        /// </summary>
        /// <param name="proteinAccession">the protein acession</param>
        /// <param name="peptides">list of all peptides</param>
        /// <param name="groups">list of all groups</param>
        /// <param name="reqNumUnmodPeptides">min num of unmodified peptides that must be observed for a protein in order to consider the protein group (default=1)</param>
        /// <param name="reqNumModPeptides">min num of modified peptides that must be observed for a protein in order to consider the protein group (default=3)</param>
        /// <param name="reqNumOfPepeptides">min num of peptides that must be observed for a protein in order to consider the protein group </param>
        /// <param name="correlationCutOff">min value at which two peptides will be considered to be correlated</param>
        /// <param name="compareUnmod">if false (default) only compare modified peptides to baseline, not unmodified peptides</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        /// <param name="groupToCompare">single group to compare all other groups to (not a required parameter)</param>
        /// <param name="useProtBaseline">correct useProt bool</param>
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
                    new Peptide("Seq1[a]", "Seq1", "Prot", "Gene", "Organism",
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

        /// <summary>
        /// Test to check that the correct number of pairwise compairisons are being calculated
        /// </summary>
        /// <param name="proteinAccession">the protein acession</param>
        /// <param name="peptides">list of all peptides</param>
        /// <param name="groups">list of all groups</param>
        /// <param name="groupToCompare">single group to compare all other groups to (not a required parameter)</param>
        /// <param name="numPairwiseCompairisonsBaselineAllGroupProt">correct number of protein pairwise compairisons when comparing all groups</param>
        /// <param name="numPairwiseCompairisonsBaselineSetGroupProt">correct number of protein pairwise compairisons when comparing against one group</param>
        /// <param name="numPairwiseCompairisonsBaselineAllGroupPTM">correct number of ptm pairwise compairisons when comparing all groups</param>
        /// <param name="numPairwiseCompairisonsBaselineSetGroupPTM">correct number of ptm pairwise compairisons when comparing against one group</param>
        [Test]
        [TestCaseSource("_calCompairison")]
        public void ProteinGroup_calcComparison_Pass(string proteinAccession, List<Peptide> peptides, List<string> groups, string groupToCompare,
            int numPairwiseCompairisonsBaselineAllGroupProt, int numPairwiseCompairisonsBaselineSetGroupProt, int numPairwiseCompairisonsBaselineAllGroupPTM,
            int numPairwiseCompairisonsBaselineSetGroupPTM)
        {

            ProteinGroup ProteinGroupBaselineAllGroupTest = new ProteinGroup(proteinAccession, peptides, groups, 3, 1, 4,
                3, 3, 0.5, false, 3);
            ProteinGroup ProteinGroupBaselineSetGroupTest = new ProteinGroup(proteinAccession, peptides, groups, 3, 1, 4,
                3, 3, 0.5, false, 3, groupToCompare);

            Assert.AreEqual(numPairwiseCompairisonsBaselineAllGroupProt, ProteinGroupBaselineAllGroupTest.ProteinPairwiseComparisons.Count());
            Assert.AreEqual(numPairwiseCompairisonsBaselineSetGroupProt, ProteinGroupBaselineSetGroupTest.ProteinPairwiseComparisons.Count());
            Assert.AreEqual(numPairwiseCompairisonsBaselineAllGroupPTM, ProteinGroupBaselineAllGroupTest.PTMPairwiseCompairisons.Count());
            Assert.AreEqual(numPairwiseCompairisonsBaselineSetGroupPTM, ProteinGroupBaselineSetGroupTest.PTMPairwiseCompairisons.Count());
        }
    }
}

