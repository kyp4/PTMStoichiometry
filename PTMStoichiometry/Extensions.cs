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

    }
}
