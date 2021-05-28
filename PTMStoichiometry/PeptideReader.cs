using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace PTMStoichiometry20210414a
{   
    //class to read in Peptide data from FlashLFQ data
    //PeptideReader class written with Dr. Shortreed
    public class PeptideReader
    {
        public static List<Peptide> ReadTsv(string peptidefilepath, string groupfilepath, int reqNumPepMeasurements)
        {
            List<Peptide> peptides = new List<Peptide>();
            string[] lines = File.ReadAllLines(peptidefilepath, Encoding.UTF8);
            string[] files = GetFile(lines[0]);
            Dictionary<string, string> groups = GetGroups(groupfilepath);
            List<string> groupList = GetGroupList(groupfilepath);

            for (int i = 1; i < lines.Length; i++)
            {
                peptides.Add(new Peptide(lines[i], files, groups, groupList, reqNumPepMeasurements));
            }
            return peptides;
        }

        public static Dictionary<string, string> GetGroups(string groupfilepath)
        {
            string[] groupsWithFiles = File.ReadAllLines(groupfilepath, Encoding.UTF8);
            Dictionary<string, string> groups = new Dictionary<string, string>();
            for (int i = 0; i < groupsWithFiles.Length; i++)
            {
                string[] temp = groupsWithFiles[i].Split('\t');
                groups.Add(temp[0], temp[1]);
            }

            return groups;
        }

        public static List<string> GetGroupList(string groupfilepath)
        {
            string[] groupsWithFiles = File.ReadAllLines(groupfilepath, Encoding.UTF8);
            List<string> groups = new List<string>();
            for (int i = 0; i < groupsWithFiles.Length; i++)
            {
                string[] temp = groupsWithFiles[i].Split('\t');
                groups.Add(temp[1]);
            }

            return groups.Distinct().ToList();
        }

        private static string[] GetFile(string v)
        {
            string[] headerEntries = v.Split('\t');
            int length = (headerEntries.Length - 5) / 2;
            for (int i = 5; i < length + 5; i++)
            {
                headerEntries[i] = headerEntries[i].Replace("Intensity_", "");
            }
            return headerEntries.SubArray(5, length);
        }
    }
}