using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace PTMStoichiometry
{   
    //class to read in Peptide data from FlashLFQ data
    //PeptideReader class written with Dr. Shortreed
    public class PeptideReader
    {
        /// <summary>
        /// Function to read in Peptide data from tsv files
        /// </summary>
        /// <param name="peptidefilepath">file path to the tsv file containing Peptide data</param>
        /// <param name="groupfilepath">file path to tsv file containing group data</param>
        /// <param name="reqNumPepMeasurements">min num of peptide intensities that must be observed 
        /// (non zero -> MS or MSMS detection)</param>
        /// <param name="intensityIndex">index where Intensity data starts in file</param>
        /// <param name="dataType">type of data: FlashLFQ or MaxQuant</param>
        /// <returns>list of peptides</returns>
        public static List<Peptide> ReadTsv(string peptidefilepath, string groupfilepath, int reqNumPepMeasurements,
            int intensityIndex, string dataType)
        {
            List<Peptide> peptides = new List<Peptide>();
            string[] lines = File.ReadAllLines(peptidefilepath, Encoding.UTF8);
   
            Dictionary<string, string> groups = GetGroups(groupfilepath);
            List<string> groupList = GetGroupList(groupfilepath);

            for (int i = 1; i < lines.Length; i++)
            {
                peptides.Add(new Peptide(lines[i], groups, groupList, reqNumPepMeasurements, intensityIndex, dataType));
            }
            return peptides;
        }

        /// <summary>
        /// Function to build dictionary of which file is related to which group 
        /// </summary>
        /// <param name="groupfilepath">file path to tsv file containing group data</param>
        /// <returns>dictionary of which file is in which group</returns>
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

        /// <summary>
        /// Function which returns a list all groups
        /// </summary>
        /// <param name="groupfilepath">file path to tsv file containing group data</param>
        /// <returns>list of all groups</returns>
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

        /// <summary>
        /// Function to find index of a key string in a header
        /// </summary>
        /// <param name="filepath">file path to the tsv file with a header containing find</param>
        /// <param name="find">string to find in header</param>
        /// <returns>the index of find in the header</returns>
        public static int IndexFind(string filepath, string find)
        {
            string headerline = File.ReadAllLines(filepath, Encoding.UTF8)[0];
            List<string> header = headerline.Split("\t").ToList();

            return header.IndexOf(find);
        }
    }
}