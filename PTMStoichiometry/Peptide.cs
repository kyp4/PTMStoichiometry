using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace PTMStoichiometry
{
    //class to work with peptide data from MM and FlashLFQ
    public class Peptide
    {
        //Seq w mods from MM search output
        public string Sequence { get; }
        //Seq without mods 
        public string BaseSeq { get; }
        //Seq mods
        public List<PostTranslationalModification> PostTranslationalModifications { get; }
        //Accession number(s)
        public List<string> ProteinGroup { get; } 
        public string GeneName { get; }
        //list of intensities
        public List<Intensity> Intensities { get; }
        // is true is peptide is in only one protein group, false otherwise
        public bool IsUnique { get; }
        // is true if peptide is detected (>0) in all group the min num of times (reqNumPepMeasurements, default=3), false otherwise
        public bool DetectedMinNum { get; }
        //true if peptide is modified and modification is not one of the "to ignore" modifications
        public bool Mod { get; } 

        //list of fixed PTMs to exclude as PTMs
        IEnumerable<string> fixedPTMs = File.ReadLines(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\PTMStoichiometry\FixedPTMs.txt");

        //list of roman numerals in [] to get rid of [] so can pull out PTMs between []
        IEnumerable<string> romanNumeralReplace = File.ReadLines(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\PTMStoichiometry\RomanNumeralReplace.txt");
        
        /// <summary>
        /// Object to store Peptide data
        /// </summary>
        /// <param name="Seq">Peptide sequence with modifications</param>
        /// <param name="BaseSeq">Peptide base sequence</param>
        /// <param name="ProteinGroup">Protein group the peptide belongs to</param>
        /// <param name="GeneName">Gene the peptide belongs to</param>
        /// <param name="Organism">Organism the peptide belongs to</param>
        /// <param name="Intensities">Intensities measured for the peptide</param>
        /// <param name="groupsList">Groups this peptide was measured for</param>
        /// <param name="reqNumPepMeasurements">number of intensity measurements required for object to be usable</param>
        public Peptide(string Seq, string BaseSeq, string ProteinGroup, string GeneName, string Organism, List<Intensity> Intensities, 
            List<string> groupsList, int reqNumPepMeasurements)
        {
            this.Sequence = Seq;
            this.BaseSeq = BaseSeq;
            this.PostTranslationalModifications = GetModifications(this.Sequence);
            this.ProteinGroup = ProteinGroup.Split(";").ToList();
            this.IsUnique = this.ProteinGroup.Count() == 1;
            this.GeneName = GeneName;
            this.Intensities = Intensities;
            this.DetectedMinNum = DetectCount(this.Intensities, groupsList, reqNumPepMeasurements);
            this.Mod = GetMod(this.Sequence, this.BaseSeq);
        }

        /// <summary>
        /// Object to store Peptide data coming from FlashLFQ and MaxQuant
        /// </summary>
        /// <param name="line">tab deliminated line read in from MaxQuant or FlashLFQ</param>
        /// <param name="groups">dictionary linking files to their group</param>
        /// <param name="groupsList">list of all groups</param>
        /// <param name="reqNumPepMeasurements">number of intensity measurements required for object to be usable</param>
        /// <param name="intensityIndex">index of where intensity values start in the line</param>
        /// <param name="dataType">FlashLFQ and MaxQuant</param>
        public Peptide(string line, Dictionary<string, string> groups, List<string> groupsList, int reqNumPepMeasurements,
            int intensityIndex, string dataType)
        {
            
            var spl = line.Split('\t');

            if (dataType == "FlashLFQ")
            {
                this.Sequence = spl[0];
                this.BaseSeq = spl[1];
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


        /// <summary>
        /// Function which returns a list of the post transliational modifications in a sequence
        /// </summary>
        /// <param name="sequence">peptide sequence</param>
        /// <returns>a list of the post transliational modifications in sequence</returns>
        private List<PostTranslationalModification> GetModifications(string sequence)
        {
            foreach (string num in romanNumeralReplace)
            {
                string[] numList = num.Split("\t");
                sequence = sequence.Replace(numList[0], numList[1]);
            }
            foreach (string ptm in fixedPTMs)
            {
                sequence = sequence.Replace(ptm, "");
            }

            List<PostTranslationalModification> ptms = new List<PostTranslationalModification>();
            List<string> mods = sequence.Split('[', ']').Where((item, index) => index % 2 != 0).ToList();
            List<string> localizedMods = Peptide.GetLocalizedModifications(sequence, mods);
            for (int i = 0; i < mods.Count(); i++)
            {
                ptms.Add(new PostTranslationalModification(mods[i], localizedMods[i], sequence));
            }
            return ptms;
        }

        /// <summary>
        /// Function which returns a list of strings which are the base sequence with only one modification
        /// at a time in the correct location for comparison across peptides
        /// </summary>
        /// <param name="seq">peotide sequence</param>
        /// <param name="mods">string list of modifications</param>
        /// <returns>list of strings which are the base sequence with only one modification
        /// at a time in the correct location for comparison across peptides</returns>
        private static List<string> GetLocalizedModifications(string seq, List<string> mods)
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



        /// <summary>
        /// Function to determine if the number of intesities detected is sufficent
        /// </summary>
        /// <param name="intensities">detected intensities</param>
        /// <param name="groupsList">list of groups</param>
        /// <param name="reqNumPepMeasurements">number of intensity measurements required for object to be usabl</param>
        /// <returns>true if all groups have at least reqNumPepMeasurements intensities above 0, false otherwise</returns>
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

        /// <summary>
        /// Function to determine if sequence is modified
        /// </summary>
        /// <param name="Sequence">peptide sequence</param>
        /// <param name="BaseSeq">base peptide sequence</param>
        /// <returns>true if the peptide is modified, false otherwise</returns>
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


        /// <summary>
        /// Function to take intensities from a line (FlashLFQ or MaxQuant) and place in intensities objects
        /// </summary>
        /// <param name="vs">line</param>
        /// <param name="groups">dictionary with file and group info</param>
        /// <returns>list of Intensity objects</returns>
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