using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Statistics;

namespace PTMStoichiometry20210414a
{
    //work with protein data from MM and FlashLFQ
    public class ProteinGroup
    {
        public List<Peptide> PeptidesInProtein { get; }
        public string ProteinName { get; }
        public Boolean useProt { get; set; } //determine whether prot is useful for stoichiometry calc - needs to have both mod & unmod state
        public int NumPeptidesInProtein { get; }
        public List<Peptide> BaselinePeptides { get; } //peptides to compare others to: must be unmodified, covary with each other, and not be in other proteins & must be the same baseline for all groups
        public List<PairwiseComparison> ProteinPairwiseComparisons { get; } //compare peptides by group within protein

        //reqNumOfPepeptides - min num of peptides that must be observed for a protein in order to consider it
        //reqNumModPeptides - min num of modified peptides that must be observed for a protein in order to consider it (default=1)
        //reqNumUnmodPeptides - min num of modified peptides that must be observed for a protein in order to consider it (default=3)
        //Baseline params:
        //useBaselinePeptides - if true (default) use an averaged baseline of covarying peptides
        //reqNumBaselinePeptides - min num of baseline peptides that must be observed for a protein in order to consider it (default=3)
        //correlationCutOff - min value at which two peptides will be considered to be correlated
        //compareUnmod - if false (default) only compare modified peptides to baseline, not unmodified peptides

        public ProteinGroup(string proteinAccession, List<Peptide> peptides, int reqNumUnmodPeptides, int reqNumModPeptides, int reqNumOfPepeptides) 
        {
            this.ProteinName = proteinAccession;
            this.PeptidesInProtein = peptides.Where(p => p.ProteinGroup == ProteinName).Where(p => p.DetectedMinNum).ToList();
            this.NumPeptidesInProtein = PeptidesInProtein.Count();
            this.useProt = isProteinUseful(peptides, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides);
        }

        public ProteinGroup(string proteinAccession, List<Peptide> peptides, List<string> groups, int reqNumUnmodPeptides, int reqNumModPeptides, int reqNumOfPepeptides,
        Boolean useBaselinePeptides, int reqNumBaselinePeptides, double correlationCutOff, Boolean compareUnmod, int minNumStoichiometries) : this (proteinAccession, peptides, reqNumUnmodPeptides, 
        reqNumModPeptides, reqNumOfPepeptides)
        {
            
            if (!this.useProt) return;
            
            //use averaged baseline
            if (useBaselinePeptides)
            {
                this.BaselinePeptides = getBaseLinePeptides(this.PeptidesInProtein, peptides, groups, correlationCutOff);
                if (this.BaselinePeptides.Count() < reqNumBaselinePeptides)
                {
                    this.setUseProt(false);
                    return;
                }
                this.ProteinPairwiseComparisons = calcComparison(this.PeptidesInProtein, groups, compareUnmod, minNumStoichiometries);
            }
            //don't use averaged baseline - compare each protein to others
            else
            {
                this.ProteinPairwiseComparisons = calcComparison(this.PeptidesInProtein, groups, minNumStoichiometries);
            }
        }

        //overload if given groupToCompare - single group to compare against, this is the group name (e.g. a control group) (default = null)
        public ProteinGroup(string proteinAccession, List<Peptide> peptides, List<string> groups, int reqNumUnmodPeptides, int reqNumModPeptides, int reqNumOfPepeptides,
        Boolean useBaselinePeptides, int reqNumBaselinePeptides, double correlationCutOff, Boolean compareUnmod, int minNumStoichiometries, string groupToCompare) : this(proteinAccession, peptides, reqNumUnmodPeptides,
        reqNumModPeptides, reqNumOfPepeptides)
        {

            if (!this.useProt) return;

            //use averaged baseline
            if (useBaselinePeptides)
            {
                this.BaselinePeptides = getBaseLinePeptides(this.PeptidesInProtein, peptides, groups, correlationCutOff);
                if (this.BaselinePeptides.Count() < reqNumBaselinePeptides)
                {
                    this.setUseProt(false);
                    return;
                }
                this.ProteinPairwiseComparisons = calcComparison(this.PeptidesInProtein, groups, compareUnmod, minNumStoichiometries, groupToCompare);
            }
            //don't use averaged baseline - compare each protein to others
            else
            {
                this.ProteinPairwiseComparisons = calcComparison(this.PeptidesInProtein, groups, minNumStoichiometries, groupToCompare);
            }
        }


        //function to find baseline peptides to compare against
        private List<Peptide> getBaseLinePeptides(List<Peptide> peptidesInProtein, List<Peptide> Allpeptides, List<string> groups, double correlationCutOff)
        {
            List<Peptide>  unmodPep = peptidesInProtein.Where(p => p.Sequence == p.BaseSeq).Where(p => p.IsUnique).ToList();

            //if pos Corr(X,Y) and pos Corr(X,Z) => pos Corr(Y,Z)
            List<Peptide> unmodPepCov = new List<Peptide>();
            for (int p1 = 0; p1 < unmodPep.Count(); p1++) 
            {
                List<Peptide> temp = new List<Peptide>();
                temp.Add(unmodPep[p1]);
                for (int p2 = p1+1; p2 < unmodPep.Count(); p2++)
                {
                    if (GroupCorrelation(unmodPep[p1].Intensities, unmodPep[p2].Intensities, groups) > correlationCutOff) 
                        temp.Add(unmodPep[p2]);
                }

                if (temp.Count() > unmodPepCov.Count())
                {
                    unmodPepCov = temp;
                }
            }

            return unmodPepCov;
        }

        //function to calc Correlation in each group - returns min correlation
        private double GroupCorrelation(List<Intensity> Peptide1, List<Intensity> Peptide2, List<string> groups)
        {
            List<Intensity> Peptide1Vals = Peptide1.Where(p => p.Detection == DetectionMS.MS || p.Detection == DetectionMS.MSMS).ToList();
            List<Intensity> Peptide2Vals = Peptide2.Where(p => p.Detection == DetectionMS.MS || p.Detection == DetectionMS.MSMS).ToList();

            foreach (string group in groups)
            {
                //require at least three measurements in each group
                if (Peptide1Vals.Where(p => p.GroupID == group).Count() < 3 || Peptide2Vals.Where(p => p.GroupID == group).Count() < 3) 
                {
                    return -3;
                }
            }

            double correlation = 3; 
            foreach (string group in groups)
            {
                IEnumerable<double> Pep1group = Peptide1Vals.Where(p => p.GroupID == group).Select(p => p.IntensityVal);
                IEnumerable<double> Pep2group =  Peptide2Vals.Where(p => p.GroupID == group).Select(p => p.IntensityVal);
                double temp = Correlation.Pearson(Pep1group, Pep2group);
                if (temp < correlation)
                {
                    correlation = temp;
                }
            }

            return correlation;
        }

        public void setUseProt(Boolean newBool)
        {
            this.useProt = newBool;
        }



        //function to check whether the protein meets criteria input
        private Boolean isProteinUseful(List<Peptide> pepsInProt, int reqNumUnmodPeptides, int reqNumModPeptides, int reqNumOfPepeptides)
        {
            //ensure have enough peptides, mod, & unmod peptides
            if (pepsInProt.Count() >= reqNumOfPepeptides || pepsInProt.Where(p => p.BaseSeq == p.Sequence).ToList().Count() >= reqNumUnmodPeptides || 
                pepsInProt.Where(p => p.BaseSeq != p.Sequence).ToList().Count() >= reqNumModPeptides)
            {
                return true;
            }
            return false;
        }

        //function to get list of intensities of baseline peptides
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

        //baseline case: function to compare all modified peptides in protein against baseline using PairwiseCompairison 
        private List<PairwiseComparison> calcComparison(List<Peptide> PepsInProt, List<string> groups, Boolean compareUnmod, int minNumStoichiometries)
        {
            List<PairwiseComparison> comparePeps = new List<PairwiseComparison>();
            List<Peptide> PeptidesToCompare = PepsInProt;
            //if not considering unmod peptides -> remove
            if (!compareUnmod)
            {
                PeptidesToCompare = PeptidesToCompare.Where(p => p.BaseSeq != p.Sequence).ToList();
            }
            //loop over all peptide, group, group permutations 
            for (int p1 = 0; p1 < PeptidesToCompare.Count(); ++p1)
            {
                 for (int g1 = 0; g1 < groups.Count(); ++g1)
                 {
                     for (int g2 = (g1 + 1); g2 < groups.Count(); ++g2)
                     {
                            PairwiseComparison temp = new PairwiseComparison(PeptidesToCompare[p1], getBaselineIntensities(), groups[g1], groups[g2], minNumStoichiometries);
                            if (temp.PeptideStoichiometriesGroupOne.Count() >= minNumStoichiometries && temp.PeptideStoichiometriesGroupTwo.Count() >= minNumStoichiometries) 
                            {
                                comparePeps.Add(temp);
                            }
                     }
                 }
            }
            return comparePeps;
        }

        //overload calcComparison (baseline averaged comparing against only one group) function to compare all modified peptides in protein against baseline using PairwiseCompairison 
        private List<PairwiseComparison> calcComparison(List<Peptide> PepsInProt, List<string> groups, Boolean compareUnmod, int minNumStoichiometries, string groupToCompare)
        {
            List<PairwiseComparison> comparePeps = new List<PairwiseComparison>();
            List<Peptide> PeptidesToCompare = PepsInProt;
            //if not considering unmod peptides -> remove
            if (!compareUnmod)
            {
                PeptidesToCompare = PeptidesToCompare.Where(p => p.BaseSeq != p.Sequence).ToList();
            }
            groups.Remove(groupToCompare);
            //loop over all peptide, group, groupToCompare permutations 
            for (int p1 = 0; p1 < PeptidesToCompare.Count(); ++p1)
            {
                for (int g1 = 0; g1 < groups.Count(); ++g1)
                {
                    PairwiseComparison temp = new PairwiseComparison(PeptidesToCompare[p1], getBaselineIntensities(), groups[g1], groupToCompare, minNumStoichiometries);
                    if (temp.PeptideStoichiometriesGroupOne.Count() >= minNumStoichiometries && temp.PeptideStoichiometriesGroupTwo.Count() >= minNumStoichiometries)
                    {
                        comparePeps.Add(temp);
                    }
                }
            }
            return comparePeps;
        }

        //overload calcComparison (this one is for the non baseline case comparing all groups) PairwiseCompairsons for each pair of peptides and every pair of groups
        private List<PairwiseComparison> calcComparison(List<Peptide> PepsInProt, List<string> groups, int minNumStoichiometries)
        {
            List<PairwiseComparison> comparePeps = new List<PairwiseComparison>();
            //loop over all peptide, peptide, group, group permutations
            for (int p1 = 0; p1 < this.PeptidesInProtein.Count(); ++p1)
            {
                for (int p2 = (p1+1); p2 < this.PeptidesInProtein.Count(); ++p2)
                {
                    for (int g1 = 0; g1 < groups.Count(); ++g1)
                    {
                        for (int g2 = (g1+1); g2 < groups.Count(); ++g2)
                        {
                            PairwiseComparison temp = new PairwiseComparison(PepsInProt[p1], PepsInProt[p2], groups[g1], groups[g2], minNumStoichiometries);
                            if (temp.PeptideStoichiometriesGroupOne.Count() >= minNumStoichiometries && temp.PeptideStoichiometriesGroupTwo.Count() >= minNumStoichiometries)
                            {
                                comparePeps.Add(temp);
                            }
                        }
                    }
                }
            }
            return comparePeps;
        }

        //overload calcComparison (this one is for the non baseline case comparing to one group) PairwiseCompairsons for each pair of peptides and each group to one group
        private List<PairwiseComparison> calcComparison(List<Peptide> PepsInProt, List<string> groups, int minNumStoichiometries, string groupToCompare)
        {
            List<PairwiseComparison> comparePeps = new List<PairwiseComparison>();
            groups.Remove(groupToCompare);
            //loop over all peptide, peptide, group, groupToCompare permutations
            for (int p1 = 0; p1 < this.PeptidesInProtein.Count(); ++p1)
            {
                for (int p2 = (p1 + 1); p2 < this.PeptidesInProtein.Count(); ++p2)
                {
                    for (int g1 = 0; g1 < groups.Count(); ++g1)
                    {
                        PairwiseComparison temp = new PairwiseComparison(PepsInProt[p1], PepsInProt[p2], groups[g1], groupToCompare, minNumStoichiometries);
                        if (temp.PeptideStoichiometriesGroupOne.Count() >= minNumStoichiometries && temp.PeptideStoichiometriesGroupTwo.Count() >= minNumStoichiometries)
                        {
                            comparePeps.Add(temp);
                        }
                    }
                }
            }
            return comparePeps;
        }

    }
}
