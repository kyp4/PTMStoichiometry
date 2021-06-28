using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace PTMStoichiometry
{
    //work with peptide data from MM and FlashLFQ
    public class Peptide
    {
        public string Sequence { get; } //Seq w mods from MM search output
        public string BaseSeq { get; } //Seq without mods 
        public List<PostTranslationalModification> PostTranslationalModifications { get; } //Seq mods
        public List<string> ProteinGroup { get; } //Accession number(s)
        public string GeneName { get; }
        public List<Intensity> Intensities { get; } 
        public bool IsUnique { get; } // is true is peptide is in only one protein group, false otherwise
        public bool DetectedMinNum { get; } // is true if peptide is detected (>0) in all group the min num of times (reqNumPepMeasurements, default=3), false otherwise
        public bool Mod { get; } //true if peptide is modified and modification is not one of the "to ignore" modifications

        IEnumerable<string> fixedPTMs = File.ReadLines(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\PTMStoichiometry\FixedPTMs.txt");
        public Peptide(string Seq, string BaseSeq, string ProteinGroup, string GeneName, string Organism, List<Intensity> Intensities, 
            List<string> groupsList, int reqNumPepMeasurements)
        {
            this.Sequence = Seq;
            this.BaseSeq = BaseSeq;
            this.ProteinGroup = ProteinGroup.Split(";").ToList();
            this.IsUnique = this.ProteinGroup.Count() == 1;
            this.GeneName = GeneName;
            this.Intensities = Intensities;
            this.DetectedMinNum = DetectCount(this.Intensities, groupsList, reqNumPepMeasurements);
            this.Mod = GetMod(this.Sequence, this.BaseSeq);
        }


        //read from FlashLFQ and MaxQuant
        public Peptide(string line, Dictionary<string, string> groups, List<string> groupsList, int reqNumPepMeasurements,
            int intensityIndex, string dataType)
        {
            
            var spl = line.Split('\t');

            if (dataType == "FlashLFQ")
            {
                this.Sequence = spl[0];
                this.BaseSeq = spl[1];
                //this.Modifications = this.Sequence.Split('[', ']').Where((item, index) => index % 2 != 0).ToList();
                //this.LocalizedMods = GetLocalizedModifications(this.Sequence, this.Modifications, "FlashLFQ"); 
                this.PostTranslationalModifications = GetModifications(this.Sequence);

                this.ProteinGroup = spl[2].Split(";").ToList();
                this.IsUnique = this.ProteinGroup.Count() == 1;
                this.GeneName = spl[3];
                this.Intensities = GetIntensities(spl.SubArray(intensityIndex, spl.Length - intensityIndex), groups);
                this.DetectedMinNum = DetectCount(this.Intensities, groupsList, reqNumPepMeasurements);
                this.Mod = GetMod(this.Sequence, this.BaseSeq);
            }
            else if (dataType == "MaxQuant")
            {
                this.Sequence = spl[0] + spl[1];
                this.BaseSeq = spl[0];
                this.PostTranslationalModifications = GetModifications(this.Sequence);
                this.ProteinGroup = spl[5].Split(";").ToList();
                this.IsUnique = this.ProteinGroup.Count() == 1;
                this.GeneName = spl[6];
                this.Intensities = GetIntensities(spl.SubArray(intensityIndex, spl.Length - intensityIndex), groups);
                this.DetectedMinNum = DetectCount(this.Intensities, groupsList, reqNumPepMeasurements);
                //this.Mod = this.Modifications[0] != "Unmodified";

            }
            else
            {
                throw new Exception("Error: Unknown data input type: must be FlashLFQ or MaxQuant");
            }
        }

        public static List<PostTranslationalModification> GetModifications(string sequence)
        {
            List<PostTranslationalModification> ptms = new List<PostTranslationalModification>();
            List<string> mods = sequence.Split('[', ']').Where((item, index) => index % 2 != 0).ToList();
            List<string> localizedMods = Peptide.GetLocalizedModifications(sequence, mods);
            for (int i = 0; i < mods.Count(); i++)
            {
                ptms.Add(new PostTranslationalModification(mods[i], localizedMods[i]));
            }
            return ptms;
        }

        public static List<string> GetLocalizedModifications(string seq, List<string> mods)
        {
            List<string> localized = new List<string>();
                List<string> splitSeq = seq.Split('[', ']').ToList();
                List<string> modsCopy = new List<string>(mods);

                // List<string> splitBaseSeq = splitSeq.Where((item, index) => index % 2 != 1).ToList();

                for (int i = 0; i < mods.Count(); i++)
                {
                    //if (!mods[i].Contains("on X"))
                    //{
                        modsCopy.Remove(mods[i]);
                        localized.Add(string.Join("", seq.Split('[', ']').ToList().Where(p => !modsCopy.Contains(p)).ToList()));

                        modsCopy.Add(mods[i]);
                        //localized.Add(seq.Split('[', ']').ToList().RemoveAll(p => p.StartsWith("[") if(p != mods[i])).ToString());
                   // }
                }
            
            
            return localized;
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
        private List<Intensity> GetIntensities(string[] vs, Dictionary<string, string> groups)
        {
            List<Intensity> intensities = new List<Intensity>();
            int length = groups.Count();
            int i = 0;
            foreach (KeyValuePair<string, string> file in groups)
            {
                string filename = file.Key;
                string group = file.Value;
                double intensity = Convert.ToDouble(vs[i]);
                ++i;
                intensities.Add(new Intensity(filename, group, intensity));
            }

            return intensities;
        }


    }
}