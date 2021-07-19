using System;
using System.Collections.Generic;
using System.Linq;

namespace PTMStoichiometry
{
    //helper functions
    public static class Extensions
    {

        /// <summary>
        /// Function to create subarray from an array
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="array">the original array</param> 
        /// <param name="offset">where to start copying from in the original array</param>
        /// <param name="length">the length of the final subarray</param> 
        /// <returns name="result">subarray of the array starting at offset</returns>
        public static T[] SubArray<T>(this T[] array, int offset, int length)
        {
            T[] result = new T[length];
            Array.Copy(array, offset, result, 0, length);
            return result;
        }

        /// <summary>
        /// Function to apply the Benjamini-Hochberg p-value correction
        /// calls setCorrectedpVal in PairwiseCompairison
        /// </summary>
        /// <param name="pairwiseComparisons">all the PairwiseCompirison objects which contain p-values</param>
        /// <param name="alpha">level of significance</param>
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
        /// <summary>
        /// Function to calculate the the Benjamini-Hochberg corrected p-value via tweo different methods:
        /// with the PairwiseCompairisons grouped by protein and the BH correction applied separately to 
        /// each protein, or with the BH correction applied to all PairwiseCompairisons
        /// </summary>
        /// <param name="proteins">list of protein groups</param>
        /// <param name="grouped">whether to apply BH correction separately to each protein (true)
        /// or with proteins together (false)</param>
        /// <param name="alpha">significance level</param>
        /// <see cref="BenjaminiHochberg"/>
        public static void CalcCorrectedPValue(List<ProteinGroup> proteins, bool grouped = false, double alpha = 0.05)
        {
            if (grouped)
            {
                //apply correction to the pairwise compairisons in each protein individually 
                foreach (ProteinGroup prot in proteins)
                {
                    BenjaminiHochberg(prot.ProteinPairwiseComparisons, alpha);
                    BenjaminiHochberg(prot.PTMPairwiseCompairisons, alpha);
                }
            }
            //apply correction to all pairwise compairsons at once, regardless of what protein they belong to
            else
            {
                List<PairwiseCompairison> protPairwiseComparisons = new List<PairwiseCompairison>();
                List<PairwiseCompairison> ptmPairwiseComparisons = new List<PairwiseCompairison>();
                foreach (ProteinGroup prot in proteins)
                {
                    foreach (PairwiseCompairison comp in prot.ProteinPairwiseComparisons)
                    {
                        protPairwiseComparisons.Add(comp);

                    }
                    foreach (PairwiseCompairison c in prot.PTMPairwiseCompairisons)
                    {
                        ptmPairwiseComparisons.Add(c);

                    }
                }
               
                
                BenjaminiHochberg(protPairwiseComparisons, alpha);
                BenjaminiHochberg(ptmPairwiseComparisons, alpha);
                
            }
        }

        

    }
}
