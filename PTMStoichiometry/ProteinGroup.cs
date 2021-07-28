using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using MathNet.Numerics.Statistics;

namespace PTMStoichiometry
{
    //class to work with protein data from MM and FlashLFQ
    public class ProteinGroup
    {
        //list of all peptides within a protein
        public List<Peptide> PeptidesInProtein { get; }
        //protein ascension 
        public string ProteinName { get; }
        //determine whether prot is useful for stoichiometry calc - needs to have both mod & unmod state
        public bool useProt { get; set; } 
        //number of peptides in a protein
        public int NumPeptidesInProtein { get; }
        //peptides to compare others to: must be unmodified (or fixed mods), covary with each other, 
        //and not be in other proteins & must have the same baseline for all groups
        public List<Peptide> BaselinePeptides { get; }
        //compare peptides by group within protein (peptide/baseline or peptide/peptide)
        public List<PairwiseCompairison> ProteinPairwiseComparisons { get; }
        //compare ptms by group within protein (peptides (with same ptm)/baseline or peptides (with same ptm)/peptide)
        public List<PairwiseCompairison> PTMPairwiseCompairisons { get; }

        /// <summary>
        /// object to hold data relating to a protein group
        /// </summary>
        /// <param name="proteinAccession">protein ascension </param>
        /// <param name="peptides">all peptides (not just those in protein group)</param>
        /// <param name="reqNumUnmodPeptides">min num of unmodified peptides that must be observed for a protein in order to consider the protein group (default=1)</param>
        /// <param name="reqNumModPeptides">min num of modified peptides that must be observed for a protein in order to consider the protein group (default=3)</param>
        /// <param name="reqNumOfPepeptides">min num of peptides that must be observed for a protein in order to consider the protein group </param>
        public ProteinGroup(string proteinAccession, List<Peptide> peptides, int reqNumUnmodPeptides, int reqNumModPeptides, int reqNumOfPepeptides) 
        {
            this.ProteinName = proteinAccession;
            this.PeptidesInProtein = peptides.Where(p => p.ProteinGroup.Contains(ProteinName)).Where(p => p.DetectedMinNum).ToList(); 
            this.NumPeptidesInProtein = this.PeptidesInProtein.Count();
            this.useProt = isProteinUseful(this.PeptidesInProtein, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides);
        }

        /// <summary>
        /// object to hold data relating to a protein group
        /// </summary>
        /// <param name="proteinAccession">protein ascension </param>
        /// <param name="peptides">all peptides (not just those in protein group)</param>
        /// <param name="groups">list of all groups</param>
        /// <param name="reqNumUnmodPeptides">min num of unmodified peptides that must be observed for a protein in order to consider the protein group (default=1)</param>
        /// <param name="reqNumModPeptides">min num of modified peptides that must be observed for a protein in order to consider the protein group (default=3)</param>
        /// <param name="reqNumOfPepeptides">min num of peptides that must be observed for a protein in order to consider the protein group </param>
        /// <param name="useBaselinePeptides">if true use baseline proteins (pep/baseline) else compare peptide to peptide (pep/pep)</param>
        /// <param name="reqNumBaselinePeptides">min num of baseline peptides that must be observed for a protein in order to consider it (default=3)</param>
        /// <param name="reqNumBaselineMeasurements">min num of intensities (non zero -> MS or MSMS detection) that must be observed for in a 
        ///baseline peptide (non zero -> MS or MSMS detection) - increasing this value will decrease the number of baseline peptides  
        ///that are not observed in samples and therefore the number of non numeric stoichiometry values found in baseline case</param>
        /// <param name="correlationCutOff">min value at which two peptides will be considered to be correlated</param>
        /// <param name="compareUnmod">if false (default) only compare modified peptides to baseline, not unmodified peptides</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        public ProteinGroup(string proteinAccession, List<Peptide> peptides, List<string> groups, int reqNumUnmodPeptides, int reqNumModPeptides, int reqNumOfPepeptides,
        bool useBaselinePeptides, int reqNumBaselinePeptides, int reqNumBaselineMeasurements, double correlationCutOff, bool compareUnmod, int minNumStoichiometries) : this (proteinAccession, peptides, reqNumUnmodPeptides, 
        reqNumModPeptides, reqNumOfPepeptides)
        {
            if (!this.useProt) return;
            
            //use averaged baseline
            if (useBaselinePeptides)
            {
                this.BaselinePeptides = getBaseLinePeptides(this.PeptidesInProtein, groups, correlationCutOff, reqNumBaselineMeasurements);
                if (this.BaselinePeptides.Count() < reqNumBaselinePeptides)
                {
                    this.setUseProt(false);
                    return;
                }
                this.ProteinPairwiseComparisons = calcComparison(this.PeptidesInProtein, groups, compareUnmod, minNumStoichiometries).Distinct().ToList();
                this.PTMPairwiseCompairisons = ptmcalcComparison(this.PeptidesInProtein, groups, minNumStoichiometries).Distinct().ToList(); //TODO
            }
            //don't use averaged baseline - compare each protein to others
            else
            {
                this.ProteinPairwiseComparisons = calcComparison(this.PeptidesInProtein, groups, minNumStoichiometries);
                this.PTMPairwiseCompairisons =  ptmcalcComparison(this.PeptidesInProtein, groups, minNumStoichiometries);
            }
        }

        /// <summary>
        /// object to hold data relating to a protein group
        /// </summary>
        /// <param name="proteinAccession">protein ascension </param>
        /// <param name="peptides">all peptides (not just those in protein group)</param>
        /// <param name="groups">list of all groups</param>
        /// <param name="reqNumUnmodPeptides">min num of unmodified peptides that must be observed for a protein in order to consider the protein group (default=1)</param>
        /// <param name="reqNumModPeptides">min num of modified peptides that must be observed for a protein in order to consider the protein group (default=3)</param>
        /// <param name="reqNumOfPepeptides">min num of peptides that must be observed for a protein in order to consider the protein group </param>
        /// <param name="useBaselinePeptides">if true use baseline proteins (pep/baseline) else compare peptide to peptide (pep/pep)</param>
        /// <param name="reqNumBaselinePeptides">min num of baseline peptides that must be observed for a protein in order to consider it (default=3)</param>
        /// <param name="reqNumBaselineMeasurements">min num of intensities (non zero -> MS or MSMS detection) that must be observed for in a 
        ///baseline peptide (non zero -> MS or MSMS detection) - increasing this value will decrease the number of baseline peptides  
        ///that are not observed in samples and therefore the number of non numeric stoichiometry values found in baseline case</param>
        /// <param name="correlationCutOff">min value at which two peptides will be considered to be correlated</param>
        /// <param name="compareUnmod">if false (default) only compare modified peptides to baseline, not unmodified peptides</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        /// <param name="groupToCompare">single group to compare all other groups to</param>
        public ProteinGroup(string proteinAccession, List<Peptide> peptides, List<string> groups, int reqNumUnmodPeptides, int reqNumModPeptides, int reqNumOfPepeptides,
        bool useBaselinePeptides, int reqNumBaselinePeptides, int reqNumBaselineMeasurements, double correlationCutOff, bool compareUnmod, int minNumStoichiometries, string groupToCompare) : 
            this(proteinAccession, peptides, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides)
        {
            if (!this.useProt) return;

            //use averaged baseline
            if (useBaselinePeptides)
            {
                this.BaselinePeptides = getBaseLinePeptides(this.PeptidesInProtein, groups, correlationCutOff, reqNumBaselineMeasurements);
                if (this.BaselinePeptides.Count() < reqNumBaselinePeptides)
                {
                    this.setUseProt(false);
                    return;
                }
                this.ProteinPairwiseComparisons = calcComparison(this.PeptidesInProtein, groups, compareUnmod, minNumStoichiometries, groupToCompare);
                this.PTMPairwiseCompairisons = ptmcalcComparison(this.PeptidesInProtein, groups, minNumStoichiometries, groupToCompare);
            }
            //don't use averaged baseline - compare each protein to others
            else
            {
                this.ProteinPairwiseComparisons = calcComparison(this.PeptidesInProtein, groups, minNumStoichiometries, groupToCompare);
                this.PTMPairwiseCompairisons = ptmcalcComparison(this.PeptidesInProtein, groups, minNumStoichiometries, groupToCompare);
            }
        }


        //function to find baseline peptides to compare against
        /// <summary>
        /// Function to find the baseline peptides
        /// </summary>
        /// <param name="peptidesInProtein">list of peptides in the protein</param>
        /// <param name="groups">all groups</param>
        /// <param name="correlationCutOff">min value at which two peptides will be considered to be correlated</param>
        /// <param name="reqNumBaselineMeasurements">min num of intensities (non zero -> MS or MSMS detection) that must be observed for in a 
        ///baseline peptide (non zero -> MS or MSMS detection) - increasing this value will decrease the number of baseline peptides  
        ///that are not observed in samples and therefore the number of non numeric stoichiometry values found in baseline case</param>
        /// <returns>list of baseline peptides</returns>
        private List<Peptide> getBaseLinePeptides(List<Peptide> peptidesInProtein, List<string> groups, 
            double correlationCutOff, int reqNumBaselineMeasurements)
        {
            List<Peptide>  unmodPep = peptidesInProtein.Where(p => !p.Mod).Where(p => p.IsUnique).ToList();

            //if pos Corr(X,Y) and pos Corr(X,Z) => pos Corr(Y,Z)
            List<Peptide> unmodPepCov = new List<Peptide>();
            for (int p1 = 0; p1 < unmodPep.Count(); p1++) 
            {
                List<Peptide> temp = new List<Peptide>();
                temp.Add(unmodPep[p1]);
                for (int p2 = p1+1; p2 < unmodPep.Count(); p2++)
                {
                    if (GroupCorrelation(unmodPep[p1].Intensities, unmodPep[p2].Intensities, groups, reqNumBaselineMeasurements) > correlationCutOff) 
                        temp.Add(unmodPep[p2]);
                }

                if (temp.Count() > unmodPepCov.Count())
                {
                    unmodPepCov = temp;
                }
            }

            return unmodPepCov;
        }

        /// <summary>
        /// Function to calculate the pearson correlation between the peptides within each group
        /// </summary>
        /// <param name="Peptide1">list of intensities from a peptide to compare</param>
        /// <param name="Peptide2">list of intensities from a peptide to comapre</param>
        /// <param name="groups">list of all groups</param>
        /// <param name="reqNumBaselineMeasurements">required number of measurements in each group to perform compairison (must be at least 3)</param>
        /// <returns>the min correlation between the peptide intensities within each group</returns>
        private double GroupCorrelation(List<Intensity> Peptide1, List<Intensity> Peptide2, List<string> groups, int reqNumBaselineMeasurements)
        {
            List<Intensity> Peptide1Vals = Peptide1.Where(p => p.IntensityVal > 0).ToList();
            List<Intensity> Peptide2Vals = Peptide2.Where(p => p.IntensityVal > 0).ToList();

            foreach (string group in groups)
            {
                //require at least reqNumBaselineMeasurements measurements in each group
                if (Peptide1Vals.Where(p => p.GroupID == group).Count() < reqNumBaselineMeasurements || 
                    Peptide2Vals.Where(p => p.GroupID == group).Count() < reqNumBaselineMeasurements) 
                {
                    return -3;
                }
            }

            double correlation = 3; 
            foreach (string group in groups)
            {
                IEnumerable<double> Pep1group = Peptide1.Where(p => p.GroupID == group).Select(p => p.IntensityVal);
                IEnumerable<double> Pep2group = Peptide2.Where(p => p.GroupID == group).Select(p => p.IntensityVal);
                if (Pep1group.Count() != Pep2group.Count())
                {
                    throw new Exception("Samples must be paired for baseline method- if not all groups are equal this will cause an issue with baseline method, select non-baseline method");
                }
                double temp = Correlation.Pearson(Pep1group, Pep2group);
                if (temp < correlation)
                {
                    correlation = temp;
                }
            }

            return correlation;
        }

        /// <summary>
        /// Function to set UseProt
        /// </summary>
        /// <param name="newBool">bool to set UseProt</param>
        public void setUseProt(bool newBool)
        {
            this.useProt = newBool;
        }

        /// <summary>
        /// Function to determine whether protein group meets criteria to be useful
        /// </summary>
        /// <param name="pepsInProt">peptides in protein group</param>
        /// <param name="reqNumUnmodPeptides">min num of unmodified peptides that must be observed for a protein in order to consider the protein group (default=1)</param>
        /// <param name="reqNumModPeptides">min num of modified peptides that must be observed for a protein in order to consider the protein group (default=3)</param>
        /// <param name="reqNumOfPepeptides">min num of peptides that must be observed for a protein in order to consider the protein group </param>
        /// <returns>bool indicating whether protein group meets input criteria for a useful protein group</returns>
        private bool isProteinUseful(List<Peptide> pepsInProt, int reqNumUnmodPeptides, int reqNumModPeptides, int reqNumOfPepeptides)
        {
            //ensure have enough peptides, mod, & unmod peptides
            if (pepsInProt.Count() >= reqNumOfPepeptides || pepsInProt.Where(p => !p.Mod).ToList().Count() >= reqNumUnmodPeptides || 
                pepsInProt.Where(p => p.Mod).ToList().Count() >= reqNumModPeptides)
            {
                return true;
            }
            return false;
        }

        /// <summary>
        /// function to get list of intensities of baseline peptides
        /// </summary>
        /// <returns>list of intenisties from the baseline peptides</returns>
        private List<Intensity> getBaselineIntensities()
        {
            List<Intensity> BaselinePepIntensity = new List<Intensity>(); //intensities baseline for group of interest
            foreach (Peptide basePep in this.BaselinePeptides)
            {
                foreach (Intensity i2 in basePep.Intensities)
                {
                    BaselinePepIntensity.Add(i2);
                }

            }

            return BaselinePepIntensity;
        }
 
        /// <summary>
        /// Function to calculate PairwiseCompairisons in the baseline case grouped by localization modification ratios 
        /// </summary>
        /// <param name="PepsInProt">list of peptides in protein group</param>
        /// <param name="groups">list of groups</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        /// <returns>list of pairwisecompairisons for ptms</returns>
        private List<PairwiseCompairison> ptmcalcComparison(List<Peptide> PepsInProt, List<string> groups, int minNumStoichiometries)
        {
            List<PairwiseCompairison> comparePeps = new List<PairwiseCompairison>();
            List<Peptide> PeptidesToCompare = PepsInProt.Where(p => p.Mod).ToList();

            //collect ptms
            Dictionary<string, string> ptms = new Dictionary<string, string> { };
            Dictionary<string, List<Peptide>> ptmDict = new Dictionary<string, List<Peptide>> { };
            foreach (Peptide pep in PeptidesToCompare)
            {
                //ptms.Concat(pep.PostTranslationalModifications.Select(p => p.ModificationInPeptideSequence));
                foreach (PostTranslationalModification mod in pep.PostTranslationalModifications)
                {
                    if (!ptms.ContainsKey(mod.Modification))
                    {
                        ptms.Add(mod.Modification, mod.ModificationInPeptideSequence);
                    }
                        
                    if (ptmDict.ContainsKey(mod.Modification))
                    {
                        ptmDict[mod.Modification].Add(pep);
                    }
                    else
                    {
                        ptmDict.Add(mod.Modification, new List<Peptide>() { pep });
                    }
                } 
            }

            //loop over all mod, group, group permutations 
            foreach (string ptm in ptms.Keys)
            {
                for (int g1 = 0; g1 < groups.Count(); ++g1)
                {
                    for (int g2 = (g1 + 1); g2 < groups.Count(); ++g2)
                    {
                        
                        PairwiseCompairison temp = new PairwiseCompairison(ptmDict[ptm], getBaselineIntensities(), groups[g1], groups[g2], minNumStoichiometries, ptm);
                        if (temp.PeptideStoichiometriesGroupOne.Count() >= minNumStoichiometries && temp.PeptideStoichiometriesGroupTwo.Count() >= minNumStoichiometries)
                        {
                            comparePeps.Add(temp);
                        }
                    }
                }
            }
            return comparePeps;
        }

        
        /// <summary>
        /// Function to calculate PairwiseCompairisons in the baseline case grouped by localization modification ratios
        /// with a group to compare against
        /// </summary>
        /// <param name="PepsInProt">list of peptides in protein group</param>
        /// <param name="groups">list of groups</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        /// <param name="groupToCompare"></param>
        /// <returns>list of pairwisecompairisons for ptms</returns>
        private List<PairwiseCompairison> ptmcalcComparison(List<Peptide> PepsInProt, List<string> groups, int minNumStoichiometries, string groupToCompare)
        {
            List<PairwiseCompairison> comparePeps = new List<PairwiseCompairison>();
            List<Peptide> PeptidesToCompare = PepsInProt.Where(p => p.Mod).ToList();

            //collect ptms
            Dictionary<string, string> ptms = new Dictionary<string, string> { };
            Dictionary<string, List<Peptide>> ptmDict = new Dictionary<string, List<Peptide>> { };
            foreach (Peptide pep in PeptidesToCompare)
            {
                //ptms.Concat(pep.PostTranslationalModifications.Select(p => p.ModificationInPeptideSequence));
                foreach (PostTranslationalModification mod in pep.PostTranslationalModifications)
                {
                    if (!ptms.ContainsKey(mod.Modification))
                    {
                        ptms.Add(mod.Modification, mod.ModificationInPeptideSequence);
                    }

                    if (ptmDict.ContainsKey(mod.Modification))
                    {
                        ptmDict[mod.Modification].Add(pep);
                    }
                    else
                    {
                        ptmDict.Add(mod.Modification, new List<Peptide>() { pep });
                    }
                }
            }


            //loop over all mod, group, group permutations 
            foreach (string ptm in ptms.Keys)
            {
                for (int g1 = 0; g1 < groups.Count(); ++g1)
                {

                    PairwiseCompairison temp = new PairwiseCompairison(ptmDict[ptm], getBaselineIntensities(), groups[g1], groupToCompare, minNumStoichiometries, ptm);
                    if (temp.PeptideStoichiometriesGroupOne.Count() >= minNumStoichiometries && temp.PeptideStoichiometriesGroupTwo.Count() >= minNumStoichiometries)
                    {
                        comparePeps.Add(temp);
                    }

                }
            }
            return comparePeps;
        }

         
        /// <summary>
        /// Function to calculate PairwiseCompairisons
        /// baseline case: function to compare all modified peptides in protein against baseline using PairwiseCompairison
        /// </summary>
        /// <param name="PepsInProt">list of peptides in protein group</param>
        /// <param name="groups">list of groups</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        /// <returns>list of pairwisecompairisons for peptides</returns>
        private List<PairwiseCompairison> calcComparison(List<Peptide> PepsInProt, List<string> groups, bool compareUnmod, int minNumStoichiometries)
        {
            List<PairwiseCompairison> comparePeps = new List<PairwiseCompairison>();
            List<Peptide> PeptidesToCompare = PepsInProt;
            //if not considering unmod peptides -> remove
            if (!compareUnmod)
            {
                PeptidesToCompare = PeptidesToCompare.Where(p => p.Mod).ToList();
            }
            //loop over all peptide, group, group permutations 
            for (int p1 = 0; p1 < PeptidesToCompare.Count(); ++p1)
            {
                 for (int g1 = 0; g1 < groups.Count(); ++g1)
                 {
                     for (int g2 = (g1 + 1); g2 < groups.Count(); ++g2)
                     {
                            PairwiseCompairison temp = new PairwiseCompairison(PeptidesToCompare[p1], getBaselineIntensities(), groups[g1], groups[g2], minNumStoichiometries);
                            if (temp.PeptideStoichiometriesGroupOne.Count() >= minNumStoichiometries && temp.PeptideStoichiometriesGroupTwo.Count() >= minNumStoichiometries) 
                            {
                                comparePeps.Add(temp);
                            }
                     }
                 }
            }
            return comparePeps;
        }

        /// <summary>
        /// Function to calculate PairwiseCompairisons
        /// baseline case: function to compare all modified peptides in protein against baseline using PairwiseCompairison
        /// when there is a given group to compare against
        /// </summary>
        /// <param name="PepsInProt">list of peptides in protein group</param>
        /// <param name="groups">list of groups</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        /// <param name="groupToCompare"></param>
        /// <returns>list of pairwisecompairisons for peptides</returns>
        private List<PairwiseCompairison> calcComparison(List<Peptide> PepsInProt, List<string> groups, bool compareUnmod, int minNumStoichiometries, string groupToCompare)
        {
            List<PairwiseCompairison> comparePeps = new List<PairwiseCompairison>();
            List<Peptide> PeptidesToCompare = PepsInProt;
            //if not considering unmod peptides -> remove
            if (!compareUnmod)
            {
                PeptidesToCompare = PeptidesToCompare.Where(p => p.Mod).ToList();
            }
            groups.Remove(groupToCompare);
            //loop over all peptide, group, groupToCompare permutations 
            for (int p1 = 0; p1 < PeptidesToCompare.Count(); ++p1)
            {
                for (int g1 = 0; g1 < groups.Count(); ++g1)
                {
                    PairwiseCompairison temp = new PairwiseCompairison(PeptidesToCompare[p1], getBaselineIntensities(), groups[g1], groupToCompare, minNumStoichiometries);
                    if (temp.PeptideStoichiometriesGroupOne.Count() >= minNumStoichiometries && temp.PeptideStoichiometriesGroupTwo.Count() >= minNumStoichiometries)
                    {
                        comparePeps.Add(temp);
                    }
                }
            }
            groups.Add(groupToCompare);
            return comparePeps;
        }

        /// <summary>
        /// Function to calculate PairwiseCompairisons
        /// peptide:peptide case: function to compare all modified peptides in protein against other peptidse using PairwiseCompairison
        /// </summary>
        /// <param name="PepsInProt">list of peptides in protein group</param>
        /// <param name="groups">list of groups</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        /// <returns>list of pairwisecompairisons for peptides</returns>
        private List<PairwiseCompairison> calcComparison(List<Peptide> PepsInProt, List<string> groups, int minNumStoichiometries)
        {
            List<PairwiseCompairison> comparePeps = new List<PairwiseCompairison>();
            //loop over all peptide, peptide, group, group permutations
            for (int p1 = 0; p1 < this.PeptidesInProtein.Count(); ++p1)
            {
                for (int p2 = (p1+1); p2 < this.PeptidesInProtein.Count(); ++p2)
                {
                    for (int g1 = 0; g1 < groups.Count(); ++g1)
                    {
                        for (int g2 = (g1+1); g2 < groups.Count(); ++g2)
                        {
                            PairwiseCompairison temp = new PairwiseCompairison(PepsInProt[p1], PepsInProt[p2], groups[g1], groups[g2], minNumStoichiometries);
                            
                            if (temp.PeptideStoichiometriesGroupOne.Count() >= minNumStoichiometries && 
                                temp.PeptideStoichiometriesGroupTwo.Count() >= minNumStoichiometries)
                            {
                                comparePeps.Add(temp);
                            }
                        }
                    }
                }
            }
            return comparePeps;
        }

        /// <summary>
        /// Function to calculate PairwiseCompairisons
        /// peptide:peptide case: function to compare all modified peptides in protein against other peptidse using PairwiseCompairison
        /// when a group to compare against is given
        /// </summary>
        /// <param name="PepsInProt">list of peptides in protein group</param>
        /// <param name="groups">list of groups</param>
        /// <param name="minNumStoichiometries">min num of stoichiometries req in both groups before run test</param>
        /// <param name="groupToCompare"></param>
        /// <returns>list of pairwisecompairisons for peptides</returns>
        private List<PairwiseCompairison> calcComparison(List<Peptide> PepsInProt, List<string> groups, int minNumStoichiometries, string groupToCompare)
        {
            List<PairwiseCompairison> comparePeps = new List<PairwiseCompairison>();
            
            groups.Remove(groupToCompare);
            //loop over all peptide, peptide, group, groupToCompare permutations
            for (int p1 = 0; p1 < this.PeptidesInProtein.Count(); ++p1)
            {
                for (int p2 = (p1 + 1); p2 < this.PeptidesInProtein.Count(); ++p2)
                {
                    for (int g1 = 0; g1 < groups.Count(); ++g1)
                    {
                        PairwiseCompairison temp = new PairwiseCompairison(PepsInProt[p1], PepsInProt[p2], groups[g1], groupToCompare, minNumStoichiometries);
                        if (temp.PeptideStoichiometriesGroupOne.Count() >= minNumStoichiometries && temp.PeptideStoichiometriesGroupTwo.Count() >= minNumStoichiometries)
                        {
                            comparePeps.Add(temp);
                        }
                    }
                }
            }
            groups.Add(groupToCompare);
            return comparePeps;
        }

    }
}
