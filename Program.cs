using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.Linq;

namespace PTMStoichiometry20210414a
{
    public class Program
    {
        static void Main(string[] args)
        {
            
            string filepathpeptides =  //replace with MM FlashLFQ output to run
            string filepathgroups = //replace with tab separated groups file to run
    
            Dictionary<string, string> groups = PeptideReader.GetGroups(filepathgroups);
            List<string> groupsList = PeptideReader.GetGroupList(filepathgroups);
            List<Peptide> testPeptide = PeptideReader.ReadTsv(filepathpeptides, filepathgroups);
            
            //group peptides by protein
            var proteins = testPeptide.Select(p => p.ProteinGroup).Distinct().ToArray();
            List<ProteinGroup> testProt = new List<ProteinGroup>();
            for (int i = 0; i < proteins.Length; i++)
            {
                testProt.Add(new ProteinGroup(proteins[i], testPeptide, groupsList));
            }
            List<ProteinGroup> ProteinsToUse = testProt.Where(p => p.useProt == "modandunmod").Where(p => p.CorrectedMWPValue > 0).ToList(); 
            //store protein output in XML file
            string outFile = //replace with desired output file
            XmlDocument proteinOutput = new XmlDocument();
            proteinOutput.LoadXml("<PTMStoichiometry>  </PTMStoichiometry>");
            foreach (ProteinGroup prot in ProteinsToUse) //make each protein a XML element
            {
                proteinOutput.DocumentElement.AppendChild(ProteinWriter.AddProtein(prot, proteinOutput));
            }

            proteinOutput.Save(outFile);
            
        }
    }
}
