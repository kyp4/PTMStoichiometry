using System;
using System.Collections.Generic;
using System.Text;

namespace PTMStoichiometry20210414a
{
    //Extensions class written by Dr. Shortreed
    //read in peptide data from MM
    public static class Extensions
    {
        public static T[] SubArray<T>(this T[] array, int offset, int length)
        {
            T[] result = new T[length];
            Array.Copy(array, offset, result, 0, length);
            return result;
        }

    }
}
