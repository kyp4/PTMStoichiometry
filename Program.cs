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
            
            string filepathpeptides = @"C:\Users\KAP\BioinformaticsII\2021-04-07-15-31-12_full_analysis\2021-04-07-15-31-12\Task2-SearchTask\AllQuantifiedPeptides.tsv"; //replace with MM FlashLFQ output to run
            string filepathgroups = @"C:\Users\KAP\BioinformaticsII\groupsPhosphoStudy.txt"; //replace with tab separated groups file to run
    
            Dictionary<string, string> groups = PeptideReader.GetGroups(filepathgroups);
            List<string> groupsList = PeptideReader.GetGroupList(filepathgroups);
            List<Peptide> testPeptide = PeptideReader.ReadTsv(filepathpeptides, filepathgroups);
            Extensions.IncludeSharedPeptides(testPeptide);


            //group peptides by protein
            var proteins = testPeptide.Select(p => p.ProteinGroup).Distinct().ToArray();
            List<ProteinGroup> testProt = new List<ProteinGroup>();
            for (int i = 0; i < proteins.Length; i++)
            {
                testProt.Add(new ProteinGroup(proteins[i], testPeptide, groupsList));
            }
            testProt = testProt.Where(p => p.ProteinPairwiseComparisons != null).ToList();
            Extensions.CalcCorrectedPValue(testProt, false);
            List<ProteinGroup> ProteinsToUse = testProt.Where(p => p.useProt == "modandunmod").ToList(); 
            //store protein output in XML file
            string outFile = @"C:\Users\KAP\BioinformaticsII\PTMStoichiometryFullPhosphoData6.xml"; //replace with desired output file
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
