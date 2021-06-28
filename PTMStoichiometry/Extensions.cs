using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
//using UsefulProteomicsDatabases;

namespace PTMStoichiometry
{
    //helper functions
    public static class Extensions
    {
        //read in data from tsv
        //written by Dr. Shortreed
        public static T[] SubArray<T>(this T[] array, int offset, int length)
        {
            T[] result = new T[length];
            Array.Copy(array, offset, result, 0, length);
            return result;
        }

        //apply Benjamini-Hochberg p-value correct - correct the p-values
        //set the corrected p-values in pairwise comparisons
        public static void BenjaminiHochberg(List<PairwiseCompairison> pairwiseComparisons, double alpha = 0.05)
        {
            List<PairwiseCompairison> sortedPairwiseComparisons = pairwiseComparisons.OrderBy(p => p.MWPVal).ToList();
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

            for (int i = m-2; i >= 0; i--)
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
        public static void CalcCorrectedPValue(List<ProteinGroup> protein, Boolean grouped = false, double alpha = 0.05)
        {
            if (grouped)
            {
                //apply correction to the pairwise compairisons in each protein individually 
                foreach (ProteinGroup prot in protein)
                {
                    BenjaminiHochberg(prot.ProteinPairwiseComparisons, alpha);
                }
            }
            //apply correction to all pairwise compairsons at once, regardless of what protein they belong to
            else
            {
                List<PairwiseCompairison> pairwiseComparisons = new List<PairwiseCompairison>();
                foreach (ProteinGroup prot in protein)
                {
                    foreach (PairwiseCompairison comp in prot.ProteinPairwiseComparisons)
                    {
                        pairwiseComparisons.Add(comp);
                    }
                }
               
                //foreach (ProteinGroup prot in protein)
                //{
                    BenjaminiHochberg(pairwiseComparisons, alpha);
                //}
            }
        }

        

    }
}
