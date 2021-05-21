using System;
using System.Collections.Generic;
using System.Text;

namespace PTMStoichiometry20210414a
{
    //Extensions class written with Dr. Shortreed
    public static class Extensions
    {
        //read in peptide data from MM
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

    }
}
