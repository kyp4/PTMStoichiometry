using NUnit.Framework;
using System;
using System.Collections.Generic;
using PTMStoichiometry20210414a;
using System.Linq;

namespace PTMStoichiometryTester20200415a
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
                            new Intensity("file", "group1", 10, DetectionMS.MS),
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 15, DetectionMS.MS),
                            new Intensity("file", "group1", 30, DetectionMS.MS),

                            new Intensity("file", "group2", 10, DetectionMS.MS),
                            new Intensity("file", "group2", 20, DetectionMS.MS),
                            new Intensity("file", "group2", 15, DetectionMS.MS),
                            new Intensity("file", "group2", 30, DetectionMS.MS) 
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq2", "Seq2", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 40, DetectionMS.MS),
                            new Intensity("file", "group1", 30, DetectionMS.MS),
                            new Intensity("file", "group1", 60, DetectionMS.MS),

                            new Intensity("file", "group2", 20, DetectionMS.MS),
                            new Intensity("file", "group2", 40, DetectionMS.MS),
                            new Intensity("file", "group2", 30, DetectionMS.MS),
                            new Intensity("file", "group2", 60, DetectionMS.MS) 
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
                            new Intensity("file", "group1", 10, DetectionMS.MS),
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 15, DetectionMS.MS),
                            new Intensity("file", "group1", 30, DetectionMS.MS),

                            new Intensity("file", "group2", 10, DetectionMS.MS),
                            new Intensity("file", "group2", 20, DetectionMS.MS),
                            new Intensity("file", "group2", 15, DetectionMS.MS),
                            new Intensity("file", "group2", 30, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq2", "Seq2", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 40, DetectionMS.MS),
                            new Intensity("file", "group1", 30, DetectionMS.MS),
                            new Intensity("file", "group1", 60, DetectionMS.MS),

                            new Intensity("file", "group2", 20, DetectionMS.MS),
                            new Intensity("file", "group2", 40, DetectionMS.MS),
                            new Intensity("file", "group2", 30, DetectionMS.MS),
                            new Intensity("file", "group2", 60, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq3", "Seq3", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 100, DetectionMS.MS),
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 50, DetectionMS.MS),
                            new Intensity("file", "group1", 600, DetectionMS.MS),

                            new Intensity("file", "group2", 80, DetectionMS.MS),
                            new Intensity("file", "group2", 600, DetectionMS.MS),
                            new Intensity("file", "group2", 40, DetectionMS.MS),
                            new Intensity("file", "group2", 3, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq4", "Seq4", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 40, DetectionMS.MS),
                            new Intensity("file", "group1", 30, DetectionMS.MS),
                            new Intensity("file", "group1", 60, DetectionMS.MS),

                            new Intensity("file", "group2", 20, DetectionMS.MS),
                            new Intensity("file", "group2", 40, DetectionMS.MS),
                            new Intensity("file", "group2", 30, DetectionMS.MS),
                            new Intensity("file", "group2", 60, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq4", "Seq4", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 40, DetectionMS.MS),
                            new Intensity("file", "group1", 30, DetectionMS.MS),
                            new Intensity("file", "group1", 60, DetectionMS.MS),

                            new Intensity("file", "group2", 20, DetectionMS.MS),
                            new Intensity("file", "group2", 40, DetectionMS.MS),
                            new Intensity("file", "group2", 30, DetectionMS.MS),
                            new Intensity("file", "group2", 60, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq5a", "Seq5", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 40, DetectionMS.MS),
                            new Intensity("file", "group1", 30, DetectionMS.MS),
                            new Intensity("file", "group1", 60, DetectionMS.MS),

                            new Intensity("file", "group2", 20, DetectionMS.MS),
                            new Intensity("file", "group2", 40, DetectionMS.MS),
                            new Intensity("file", "group2", 30, DetectionMS.MS),
                            new Intensity("file", "group2", 60, DetectionMS.MS)
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
                            new Intensity("file", "group1", 10, DetectionMS.MS),
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 15, DetectionMS.MS),
                            new Intensity("file", "group1", 30, DetectionMS.MS),

                            new Intensity("file", "group2", 10, DetectionMS.MS),
                            new Intensity("file", "group2", 20, DetectionMS.MS),
                            new Intensity("file", "group2", 15, DetectionMS.MS),
                            new Intensity("file", "group2", 30, DetectionMS.MS),

                            new Intensity("file", "group3", 10, DetectionMS.MS),
                            new Intensity("file", "group3", 20, DetectionMS.MS),
                            new Intensity("file", "group3", 15, DetectionMS.MS),
                            new Intensity("file", "group3", 30, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2", "group3" }, 1),
                   new Peptide("Seq2", "Seq2", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 40, DetectionMS.MS),
                            new Intensity("file", "group1", 30, DetectionMS.MS),
                            new Intensity("file", "group1", 60, DetectionMS.MS),

                            new Intensity("file", "group2", 20, DetectionMS.MS),
                            new Intensity("file", "group2", 40, DetectionMS.MS),
                            new Intensity("file", "group2", 30, DetectionMS.MS),
                            new Intensity("file", "group2", 60, DetectionMS.MS),

                            new Intensity("file", "group3", 20, DetectionMS.MS),
                            new Intensity("file", "group3", 40, DetectionMS.MS),
                            new Intensity("file", "group3", 30, DetectionMS.MS),
                            new Intensity("file", "group3", 60, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2" , "group3"}, 1),
                   new Peptide("Seq3", "Seq3", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 100, DetectionMS.MS),
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 50, DetectionMS.MS),
                            new Intensity("file", "group1", 600, DetectionMS.MS),

                            new Intensity("file", "group2", 80, DetectionMS.MS),
                            new Intensity("file", "group2", 600, DetectionMS.MS),
                            new Intensity("file", "group2", 40, DetectionMS.MS),
                            new Intensity("file", "group2", 3, DetectionMS.MS),

                            new Intensity("file", "group3", 100, DetectionMS.MS),
                            new Intensity("file", "group3", 20, DetectionMS.MS),
                            new Intensity("file", "group3", 50, DetectionMS.MS),
                            new Intensity("file", "group3", 600, DetectionMS.MS)

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
                            new Intensity("file", "group1", 10, DetectionMS.MS),
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 15, DetectionMS.MS),
                            new Intensity("file", "group1", 30, DetectionMS.MS),

                            new Intensity("file", "group2", 10, DetectionMS.MS),
                            new Intensity("file", "group2", 20, DetectionMS.MS),
                            new Intensity("file", "group2", 15, DetectionMS.MS),
                            new Intensity("file", "group2", 30, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq2", "Seq2", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 40, DetectionMS.MS),
                            new Intensity("file", "group1", 30, DetectionMS.MS),
                            new Intensity("file", "group1", 60, DetectionMS.MS),

                            new Intensity("file", "group2", 20, DetectionMS.MS),
                            new Intensity("file", "group2", 40, DetectionMS.MS),
                            new Intensity("file", "group2", 30, DetectionMS.MS),
                            new Intensity("file", "group2", 60, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq3", "Seq3", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 100, DetectionMS.MS),
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 50, DetectionMS.MS),
                            new Intensity("file", "group1", 600, DetectionMS.MS),

                            new Intensity("file", "group2", 80, DetectionMS.MS),
                            new Intensity("file", "group2", 600, DetectionMS.MS),
                            new Intensity("file", "group2", 40, DetectionMS.MS),
                            new Intensity("file", "group2", 3, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2"}, 1),
                   new Peptide("Seq4", "Seq4", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 150, DetectionMS.MS),
                            new Intensity("file", "group1", 50, DetectionMS.MS),
                            new Intensity("file", "group1", 75, DetectionMS.MS),
                            new Intensity("file", "group1", 200, DetectionMS.MS),

                            new Intensity("file", "group2", 90, DetectionMS.MS),
                            new Intensity("file", "group2", 100, DetectionMS.MS),
                            new Intensity("file", "group2", 50, DetectionMS.MS),
                            new Intensity("file", "group2", 10, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2"}, 1),
                   new Peptide("Seq5", "Seq5", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 160, DetectionMS.MS),
                            new Intensity("file", "group1", 55, DetectionMS.MS),
                            new Intensity("file", "group1", 80, DetectionMS.MS),
                            new Intensity("file", "group1", 240, DetectionMS.MS),

                            new Intensity("file", "group2", 110, DetectionMS.MS),
                            new Intensity("file", "group2", 120, DetectionMS.MS),
                            new Intensity("file", "group2", 75, DetectionMS.MS),
                            new Intensity("file", "group2", 24, DetectionMS.MS)
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
                            new Intensity("file", "group1", 10, DetectionMS.MS),
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 15, DetectionMS.MS),
                            new Intensity("file", "group1", 30, DetectionMS.MS),

                            new Intensity("file", "group2", 10, DetectionMS.MS),
                            new Intensity("file", "group2", 20, DetectionMS.MS),
                            new Intensity("file", "group2", 15, DetectionMS.MS),
                            new Intensity("file", "group2", 30, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq2", "Seq2", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 40, DetectionMS.MS),
                            new Intensity("file", "group1", 30, DetectionMS.MS),
                            new Intensity("file", "group1", 60, DetectionMS.MS),

                            new Intensity("file", "group2", 20, DetectionMS.MS),
                            new Intensity("file", "group2", 40, DetectionMS.MS),
                            new Intensity("file", "group2", 30, DetectionMS.MS),
                            new Intensity("file", "group2", 60, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq3", "Seq3", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 100, DetectionMS.MS),
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 50, DetectionMS.MS),
                            new Intensity("file", "group1", 600, DetectionMS.MS),

                            new Intensity("file", "group2", 80, DetectionMS.MS),
                            new Intensity("file", "group2", 600, DetectionMS.MS),
                            new Intensity("file", "group2", 40, DetectionMS.MS),
                            new Intensity("file", "group2", 3, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2"}, 1),
                   new Peptide("Seq4", "Seq4", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 150, DetectionMS.MS),
                            new Intensity("file", "group1", 50, DetectionMS.MS),
                            new Intensity("file", "group1", 75, DetectionMS.MS),
                            new Intensity("file", "group1", 200, DetectionMS.MS),

                            new Intensity("file", "group2", 90, DetectionMS.MS),
                            new Intensity("file", "group2", 100, DetectionMS.MS),
                            new Intensity("file", "group2", 50, DetectionMS.MS),
                            new Intensity("file", "group2", 10, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2"}, 1),
                   new Peptide("Seq5", "Seq5", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 160, DetectionMS.MS),
                            new Intensity("file", "group1", 55, DetectionMS.MS),
                            new Intensity("file", "group1", 80, DetectionMS.MS),
                            new Intensity("file", "group1", 240, DetectionMS.MS),

                            new Intensity("file", "group2", 110, DetectionMS.MS),
                            new Intensity("file", "group2", 120, DetectionMS.MS),
                            new Intensity("file", "group2", 75, DetectionMS.MS),
                            new Intensity("file", "group2", 24, DetectionMS.MS)
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
            Boolean useBaselinePeptides, int reqNumBaselinePeptides, double correlationCutOff, Boolean compareUnmod, int minNumStoichiometries, string groupToCompare, List<string> baselinePeptideSeq)
        {
            Extensions.IncludeSharedPeptides(peptides, true); //set isUnique
            ProteinGroup ProteinGroupAllGroupsTest = new ProteinGroup(proteinAccession, peptides, groups, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides, useBaselinePeptides, 
                reqNumBaselinePeptides, correlationCutOff, compareUnmod, minNumStoichiometries);
            ProteinGroup ProteinGroupSetGroupTest = new ProteinGroup(proteinAccession, peptides, groups, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides, useBaselinePeptides,
                reqNumBaselinePeptides, correlationCutOff, compareUnmod, minNumStoichiometries, groupToCompare);

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
                            new Intensity("file", "group1", 10, DetectionMS.MS),
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 15, DetectionMS.MS),
                            new Intensity("file", "group1", 30, DetectionMS.MS),

                            new Intensity("file", "group2", 10, DetectionMS.MS),
                            new Intensity("file", "group2", 20, DetectionMS.MS),
                            new Intensity("file", "group2", 15, DetectionMS.MS),
                            new Intensity("file", "group2", 30, DetectionMS.MS)
                        },
                        new List<string>() { "group1", "group2" }, 1),
                   new Peptide("Seq2", "Seq2", "Prot", "Gene", "Organism",
                        new List<Intensity>() {
                            new Intensity("file", "group1", 20, DetectionMS.MS),
                            new Intensity("file", "group1", 40, DetectionMS.MS),
                            new Intensity("file", "group1", 30, DetectionMS.MS),
                            new Intensity("file", "group1", 60, DetectionMS.MS),

                            new Intensity("file", "group2", 20, DetectionMS.MS),
                            new Intensity("file", "group2", 40, DetectionMS.MS),
                            new Intensity("file", "group2", 30, DetectionMS.MS),
                            new Intensity("file", "group2", 60, DetectionMS.MS)
                        },
                    new List<string>() { "group1", "group2" }, 1)
                },
                new List<string> { "group1", "group2" },
                3, 1, 4, 0.5, false, 3, "group1",
                false, false
            }
           
            
            
           
        };

        [Test]
        [TestCaseSource("_useProt")]
        public void ProteinGroup_useProt_Pass(string proteinAccession, List<Peptide> peptides, List<string> groups, int reqNumUnmodPeptides, int reqNumModPeptides, int reqNumOfPepeptides,
            double correlationCutOff, Boolean compareUnmod, int minNumStoichiometries, string groupToCompare, bool useProtBaseline, bool useProtPeptidePeptide)
        {
            ProteinGroup ProteinGroupBaselineAllGroupTest = new ProteinGroup(proteinAccession, peptides, groups, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides, true,
                reqNumUnmodPeptides, correlationCutOff, compareUnmod, minNumStoichiometries);
            ProteinGroup ProteinGroupBaselineSetGroupTest = new ProteinGroup(proteinAccession, peptides, groups, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides, true,
                reqNumUnmodPeptides, correlationCutOff, compareUnmod, minNumStoichiometries, groupToCompare);
            ProteinGroup ProteinGroupPeptidePeptideAllGroupTest = new ProteinGroup(proteinAccession, peptides, groups, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides, false,
                reqNumUnmodPeptides, correlationCutOff, compareUnmod, minNumStoichiometries);
            ProteinGroup ProteinGroupPeptidePeptideSetGroupTest = new ProteinGroup(proteinAccession, peptides, groups, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides, false,
                reqNumUnmodPeptides, correlationCutOff, compareUnmod, minNumStoichiometries, groupToCompare);

            Assert.AreEqual(useProtBaseline, ProteinGroupBaselineAllGroupTest.useProt);
            Assert.AreEqual(useProtBaseline, ProteinGroupBaselineSetGroupTest.useProt);
            Assert.AreEqual(useProtPeptidePeptide, ProteinGroupPeptidePeptideAllGroupTest.useProt);
            Assert.AreEqual(useProtPeptidePeptide, ProteinGroupPeptidePeptideSetGroupTest.useProt);
        }

    }
}
