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
        public string FileName { get; }
        public double IntensityVal { get; }
        public string GroupID { get; } //read from a separate tsv file with both the filenames (no Intensity_ or .raw) and their group

        public Intensity(string file, string group, double intensity)
        {
            this.FileName = file;
            this.GroupID = group;
            this.IntensityVal = intensity;
        }
    }
}
