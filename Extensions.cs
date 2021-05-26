using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;

namespace PTMStoichiometry20210414a
{
    public static class Extensions
    {
        //read in peptide data from MM
        //written by Dr. Shortreed
        public static T[] SubArray<T>(this T[] array, int offset, int length)
        {
            T[] result = new T[length];
            Array.Copy(array, offset, result, 0, length);
            return result;
        }

        //remove shared peptides if parameter set to false, otherwise include shared peptides in calculations
        //sets the IsUnique parameter in Peptide
        public static List<Peptide> IncludeSharedPeptides(List<Peptide> Peps, Boolean IncludeRazorPep = false)
        {

            foreach (Peptide pep in Peps)
            {
                pep.setIsUnique(Peps);
            }

            List<Peptide> UniquePeps = new List<Peptide>();
            if (!IncludeRazorPep)
            {
                foreach (Peptide pep in Peps)
                {
                    if (pep.IsUnique)
                    {
                        UniquePeps.Add(pep);
                    }
                }
            }

            return UniquePeps;
        }

        /*
        //apply Benjamini-Hochberg correction to p-values - correct significance
        public static double BenjaminiHochberg(List<double> pvals, double alpha = 0.5)
        {
            pvals.Sort();
            int m = pvals.Count();
            double corrected = 0;
            for (int i = 0; i < pvals.Count(); i++)
            {
                if (pvals[i] <= (i + 1) * alpha / m)
                {
                    corrected = pvals[i];
                }
                else
                {
                    return corrected;
                }
            }
            return corrected;
        }
        */

        //apply Benjamini-Hochberg p-value correct - correct the p-values
        //set the corrected p-values in pairwise comparisons
        public static void BenjaminiHochberg(List<PairwiseComparison> pairwiseComparisons, double alpha = 0.5)
        {
            List<PairwiseComparison> sortedPairwiseComparisons = pairwiseComparisons.OrderBy(p => p.MWPVal).Reverse().ToList();
            int m = sortedPairwiseComparisons.Count();
            if (m == 0)
            {
                return;
            }

            //set top p-value to set correct upper limit
            sortedPairwiseComparisons[m-1].setCorrectedpVal(sortedPairwiseComparisons[m-1].MWPVal);
            if (sortedPairwiseComparisons[m - 1].CorrectedpVal > 1)
            {
                sortedPairwiseComparisons[m - 1].setCorrectedpVal(1);
            }

            for (int i = m-2; i > 0; i--)
            {
                //pval*#pvals/#in sorted order
                sortedPairwiseComparisons[i].setCorrectedpVal(sortedPairwiseComparisons[i].MWPVal * m / (i + 1));
                //p-values must be in ascending order 
                if (sortedPairwiseComparisons[i].CorrectedpVal > sortedPairwiseComparisons[i + 1].CorrectedpVal)
                {
                    sortedPairwiseComparisons[i].setCorrectedpVal(sortedPairwiseComparisons[i + 1].CorrectedpVal);
                }
            }
        }

        //two different methods of considering peptides for p-value correction: 
        //either correction within each protein (grouped) or across all proteins
        public static void CalcCorrectedPValue(List<ProteinGroup> protein, Boolean grouped = false)
        {
            if (grouped)
            {
                //apply correction to the pairwise compairisons in each protein individually 
                foreach (ProteinGroup prot in protein)
                {
                    BenjaminiHochberg(prot.ProteinPairwiseComparisons);
                }
            }
            //apply correction to all pairwise compairsons at once, regardless of what protein they belong to
            else
            {
                List<PairwiseComparison> pairwiseComparisons = new List<PairwiseComparison>();
                foreach (ProteinGroup prot in protein)
                {
                    foreach (PairwiseComparison comp in prot.ProteinPairwiseComparisons)
                    {
                        pairwiseComparisons.Add(comp);
                    }
                }
               
                foreach (ProteinGroup prot in protein)
                {
                    BenjaminiHochberg(prot.ProteinPairwiseComparisons);
                }
            }
        }

        /*
        //two different methods of considering peptides for p-value correction: 
        //either correction within each protein (grouped) or across all proteins
        public static void CalcCorrectedPValue(List<ProteinGroup> protein, Boolean grouped = false)
        {
            if (grouped)
            {
                //apply correction to the pairwise compairisons in each protein individually 
                foreach (ProteinGroup prot in protein) 
                {
                    List<double> pvals = new List<double>();
                    foreach (PairwiseComparison comp in prot.ProteinPairwiseComparisons)
                    {
                        pvals.Add(comp.MWPVal);
                    }
                    double pval = BenjaminiHochberg(pvals);
                    prot.setCorrectedMWPVal(pval);
                }
            }
            //apply correction to all pairwise compairsons at once, regardless of what protein they belong to
            else
            {
                List<double> pvals = new List<double>();
                foreach (ProteinGroup prot in protein)
                {
                    foreach (PairwiseComparison comp in prot.ProteinPairwiseComparisons)
                    {
                        pvals.Add(comp.MWPVal);
                    }
                }
                double pval = BenjaminiHochberg(pvals);
                foreach (ProteinGroup prot in protein)
                {
                    prot.setCorrectedMWPVal(pval);
                }
            }
        }
        */

    }
}
