using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PTMStoichiometry20210414a
{
    //work with protein data from MM and FlashLFQ
    public class ProteinGroup
    {
        public List<Peptide> PeptidesInProtein { get; }
        public string ProteinName { get; }
        public string useProt { get; } //determine whether prot is useful for stochiometry calc - needs to have both mod & unmod state

        public int NumPeptidesInProtein { get; }

        public List<PairwiseComparison> ProteinPairwiseComparisons { get; } //compare peptides by group within protein
        public double CorrectedMWPValue { get; } //Benjamini-Hochberg correct p-value
        public ProteinGroup(string proteinAccession, List<Peptide> peptides, List<string> groups)
        {
            this.ProteinName = proteinAccession;
            this.PeptidesInProtein = peptides.Where(p => p.ProteinGroup == ProteinName).ToList(); //need to filter peptides - only want peptides which are in all groups and observed in 2/3 observations
            this.NumPeptidesInProtein = PeptidesInProtein.Count();
            this.useProt = isProteinUseful(PeptidesInProtein);
            if ((this.useProt == "modandunmod")) //only calc stoichiometries if there are multiple peptides both mod & unmod - TODO: consider expanding
            {
                this.ProteinPairwiseComparisons = calcComparison(groups);
                this.CorrectedMWPValue = BenjaminiHochberg(this.ProteinPairwiseComparisons);
            }
        }

        //function check whether is useful protein
        private string isProteinUseful(List<Peptide> pepsInProt)
        {
            if (pepsInProt.Count() < 2)
            {
                return "unique";
            }
            else if (pepsInProt.Where(p => p.BaseSeq == p.Sequence).ToList().Count() == 0)
            {
                return "modOnly";
            }
            else if (pepsInProt.Where(p => p.BaseSeq != p.Sequence).ToList().Count() == 0)
            {
                return "unmodOnly";
            }

            return "modandunmod";
        }

        

        //calc PairwiseCompairsons for each pair pf peptides and every pair of groups
        private List<PairwiseComparison> calcComparison(List<string> groups)
        {
            List<PairwiseComparison> comparePeps = new List<PairwiseComparison>();
            
            for (int p1 = 0; p1 < this.PeptidesInProtein.Count(); ++p1)
            {
                for (int p2 = (p1+1); p2 < this.PeptidesInProtein.Count(); ++p2)
                {
                    for (int g1 = 0; g1 < groups.Count(); ++g1)
                    {
                        for (int g2 = (g1+1); g2 < groups.Count(); ++g2)
                        {
                            PairwiseComparison temp = new PairwiseComparison(this.PeptidesInProtein[p1], this.PeptidesInProtein[p2], groups[g1], groups[g2]);
                            if (temp.PeptideStoichiometriesGroupOne.Count > 3 && temp.PeptideStoichiometriesGroupTwo.Count > 3) //p-value set to -1 if both stoichiometry groups not larger than 3 values 
                            {
                                comparePeps.Add(new PairwiseComparison(this.PeptidesInProtein[p1], this.PeptidesInProtein[p2], groups[g1], groups[g2]));
                            }
                        }
                    }
                }
            }
            return comparePeps;
        }
        public double BenjaminiHochberg(List<PairwiseComparison> Comparisons, double alpha = 0.5)
        {
            List<double> pvals = new List<double>();
            foreach (PairwiseComparison comp in Comparisons)
            {
                pvals.Add(comp.MWPVal);
            }
            //select largest p-value such that: (k+1)*alpha/#pvals > p-val, wherek = # in an ordered list of the p-vals
            pvals.Sort();
            for (int k = 0; k < pvals.Count(); ++k)
            {
                if (pvals[k] > ((k+1) * alpha) / pvals.Count())
                {
                    if (k > 0)
                    {
                        return pvals[k - 1];
                    }
                    else
                    {
                        return pvals[k];
                    }               
                }
            }
            return 0; //this makes nothing sig pvals[pvals.Count() - 1];
        }

        /*
        //calc PairwiseCompairsons for each pair pf peptides and every pair of groups
        private List<PairwiseComparison> calcComparison(Dictionary<string, string> groups)
        {
            List<PairwiseComparison> comparePeps = new List<PairwiseComparison>();
            foreach (Peptide pep1 in this.PeptidesInProtein)
            {
                foreach (Peptide pep2 in this.PeptidesInProtein)
                {
                    foreach (string g in groups.Keys)
                    {
                        foreach (string h in groups.Keys)
                        {
                            comparePeps.Add(new PairwiseComparison(pep1, pep2, groups[g], groups[h]));
                        }
                    }
                }
            }
            return comparePeps;
        }
        */
    }
}
