using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PTMStoichiometry
{
    //Works with intensities from FlashLFQ output from MM
    public class Intensity
    {
        //the file the intensity came from 
        public string FileName { get; }
        //the quantitative intensity value
        public double IntensityVal { get; }
        //read from a separate tsv file with both the filenames (no Intensity_ or .raw) and their group
        public string GroupID { get; } 

        /// <summary>
        /// Object to store intensity data
        /// </summary>
        /// <param name="file">the file the intensity data came from</param>
        /// <param name="group">the group the intensity data came from</param>
        /// <param name="intensity">the quantitative intensity value</param>
        public Intensity(string file, string group, double intensity)
        {
            this.FileName = file;
            this.GroupID = group;
            this.IntensityVal = intensity;
        }
    }
}
