using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PTMStoichiometry20210414a
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

        public Peptide(string line, string[] fileNames, Dictionary<string, string> groups)
        {
            var spl = line.Split('\t');

            this.Sequence = spl[0];
            this.BaseSeq = spl[1];
            this.ProteinGroup = spl[2];
            this.GeneName = spl[3];
            this.Organism = spl[4];
            this.Intensities = GetIntensities(spl.SubArray(5, spl.Length - 5), fileNames, groups);
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
            int count = 0;
            foreach (Peptide pep in Peps)
            {
                if (pep.Sequence == this.Sequence)
                {
                    count++;
                    if (count > 1)
                    {
                        this.IsUnique = false;
                    }
                }
            }

            this.IsUnique = true;
        }

    }
}