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
        public DetectionMS Detection { get; } //Enum regarding how well was detected and how was detected as assigned by FlashLFQ
        public string GroupID { get; } //read from a separate tsv file with both the filenames (no Intensity_ or .raw) and their group

        public Intensity(string file, string group, double intensity, DetectionMS detection)
        {
            this.FileName = file;
            this.GroupID = group;
            this.IntensityVal = intensity;
            this.Detection = detection;
        }
    }

    //FlashLFQ designations regarding peptide detection
    public enum DetectionMS { MS, MSMS, NotDetected, MSMSIdentifiedButNotQuantified, MSMSAmbiguousPeakfinding }
}
