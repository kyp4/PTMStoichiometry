using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PTMStoichiometry
{
    //work with peptide data from MM and FlashLFQ
    public class Peptide
    {
        public string Sequence { get; } //Seq w mods from MM search output
        public string BaseSeq { get; } //Seq without mods from MM search output
        public string ProteinGroup { get; } //Accession number from MM search output
        public string GeneName { get; }
        public string Organism { get; }
        public List<Intensity> Intensities { get; } 
        public Boolean IsUnique { get; set; } // is true is peptide is in only one protein group, false otherwise
        public Boolean DetectedMinNum { get; } // is true if peptide is detected (>0) in all group the min num of times (reqNumPepMeasurements, default=3), false otherwise
        public Boolean Mod { get; } //true if peptide is modified and modification is not one of the "to ignore" modifications

        IEnumerable<string> fixedPTMs = File.ReadLines(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\PTMStoichiometry\FixedPTMs.txt");
        public Peptide(string line, string[] fileNames, Dictionary<string, string> groups, List<string> groupsList, int reqNumPepMeasurements)
        {
            var spl = line.Split('\t');

            this.Sequence = spl[0];
            this.BaseSeq = spl[1];
            this.ProteinGroup = spl[2];
            this.GeneName = spl[3];
            this.Organism = spl[4];
            this.Intensities = GetIntensities(spl.SubArray(5, spl.Length - 5), fileNames, groups);
            this.DetectedMinNum = DetectCount(this.Intensities, groupsList, reqNumPepMeasurements);
            this.Mod = GetMod(this.Sequence, this.BaseSeq);
        }

        //overload for testing 
        public Peptide(string Seq, string BaseSeq, string ProteinGroup, string GeneName, string Organism, List<Intensity> Intensities, List<string> groupsList, int reqNumPepMeasurements)
        {
            this.Sequence = Seq;
            this.BaseSeq = BaseSeq;
            this.ProteinGroup = ProteinGroup;
            this.GeneName = GeneName;
            this.Organism = Organism;
            this.Intensities = Intensities;
            this.DetectedMinNum = DetectCount(this.Intensities, groupsList, reqNumPepMeasurements);
            this.Mod = GetMod(this.Sequence, this.BaseSeq);
        }

        //determine if the number of detections is greater than the min set
        private bool DetectCount(List<Intensity> intensities, List<string> groupsList, int reqNumPepMeasurements)
        {
            foreach (string group in groupsList)
            {
                if (intensities.Where(p => p.GroupID == group).Where(p => p.IntensityVal > 0).Count() < reqNumPepMeasurements)
                {
                    return false;
                }
            }
            return true;
        }

        //set mod: true if seq is mod and mod not one of the fixed mods to ignore
        private bool GetMod(string Sequence, string BaseSeq)
        {
            string seq = Sequence;
            if (Sequence == BaseSeq)
            {
                return false;
            }
            
            foreach (string ptm in fixedPTMs)
            {
                seq = seq.Replace(ptm, "");
            }
            
            if (seq == BaseSeq)
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        //take intensities from a line and place in intensities object
        private List<Intensity> GetIntensities(string[] vs, string[] fileNames, Dictionary<string, string> groups)
        {
            List<Intensity> intensities = new List<Intensity>();
            int length = vs.Count() / 2;
            for (int i = 0; i < length; i++)
            {
                string filename = fileNames[i];
                string group = groups[fileNames[i]];
                double intensity = Convert.ToDouble(vs[i]);
                DetectionMS ms = (DetectionMS)Enum.Parse(typeof(DetectionMS), vs[i + length].ToString());
                intensities.Add(new Intensity(filename, group, intensity, ms));
            }
            return intensities;
        }

        public void setIsUnique(List<Peptide> Peps)
        {
            this.IsUnique = true;
            foreach (Peptide pep in Peps)
            {
                if (pep.Sequence == this.Sequence && pep != this)
                {
                    this.IsUnique = false;
                    return;
                }
            }
            this.IsUnique = true;
        }

    }
}