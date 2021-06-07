using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.Linq;

namespace PTMStoichiometry
{
    public class Program
    {
        static void Main(string[] args)
        {
            //user parameters - will need to validate input

            //groupToCompare - single group to compare against, this is the group name (e.g. a control group) (default = null)
            //useBaselinePeptides - if true (default) use an averaged baseline of covarying peptides
            //reqNumBaselinePeptides - min num of baseline peptides that must be observed for a protein in order to consider it (default=3)
            //correlationCutOff - min value at which two peptides will be considered to be correlated
            //reqNumOfPepeptides - min num of peptides that must be observed for a protein in order to consider it
            //reqNumModPeptides - min num of modified peptides that must be observed for a protein in order to consider it (default=1)
            //reqNumUnmodPeptides - min num of modified peptides that must be observed for a protein in order to consider it (default=3)
            //reqNumPepMeasurements - min num of peptide intensities that must be observed (non zero -> MS or MSMS detection)
            //reqNumBaselineMeasurements - min num of intensities (non zero -> MS or MSMS detection) that must be observed for in a 
            //baseline peptide (non zero -> MS or MSMS detection) - increasing this value will decrease the number of baseline peptides  
            //that are not observed in samples and therefore the number of non numeric stoichiometry values found in baseline case
            //compareUnmod - if false (default) only compare modified peptides to baseline, not unmodified peptides

            //useRazorPeptides - if false (default) peptides in more than one protein are removed
            //test - currently only MW available, want to add t-test and ANOVA
            //pvalueAdjust - if true (default) use Benjamini-Hochberg analysis to adjust all p-values, else report raw p-values
            //alpha - significance value
            //minNumStoichiometries - min num of stoichiometries req in both groups before run test
            //groupPepsForPValCalc - choose to apply p-valuse correction within each protein (grouped) or across all proteins
            //alpha - chosen significance (default=0.05)

            int reqNumUnmodPeptides = 1;
            int reqNumModPeptides = 1;
            int reqNumOfPepeptides = reqNumUnmodPeptides + reqNumModPeptides;
            Boolean useBaselinePeptides = false;
            int reqNumBaselinePeptides = reqNumUnmodPeptides;
            int reqNumBaselineMeasurements = 3; //allow one missing value
            double correlationCutOff = 0.5;
            Boolean compareUnmod = false;
            int minNumStoichiometries = 3;
            int reqNumPepMeasurements = 3;
            Boolean groupPepsForPValCalc = false;
            double alpha = 0.05;



            string filepathpeptides = @"D:\PTMStoichiometry\TestData\EcoliSpikeIn\2021-06-03-11-06-06\Task3-SearchTask\AllQuantifiedPeptidesEvenGroups.txt";
            string filepathgroups = @"D:\PTMStoichiometry\TestData\EcoliSpikeIn\2021-06-03-11-06-06\Task3-SearchTask\EcoliSpikeInEvenGroups.txt"; //replace with tab separated groups file to run

            //string filepathpeptides = @"C:\Users\KAP\BioinformaticsII\2021-04-07-15-31-12_full_analysis\2021-04-07-15-31-12\Task2-SearchTask\AllQuantifiedPeptides.tsv"; //replace with MM FlashLFQ output to run
            //string filepathgroups = @"C:\Users\KAP\BioinformaticsII\groupsPhosphoStudy.txt"; //replace with tab separated groups file to run

            Dictionary<string, string> groups = PeptideReader.GetGroups(filepathgroups);
            List<string> groupsList = PeptideReader.GetGroupList(filepathgroups);
            List<Peptide> testPeptide = PeptideReader.ReadTsv(filepathpeptides, filepathgroups, reqNumPepMeasurements);
            testPeptide = Extensions.IncludeSharedPeptides(testPeptide, false);


            //group peptides by protein
            var proteins = testPeptide.Select(p => p.ProteinGroup).Distinct().ToArray();
            List<ProteinGroup> testProt = new List<ProteinGroup>();
            for (int i = 0; i < proteins.Length; i++)
            {
                testProt.Add(new ProteinGroup(proteins[i], testPeptide, groupsList, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides,
                    useBaselinePeptides, reqNumBaselinePeptides, reqNumBaselineMeasurements, correlationCutOff, compareUnmod, minNumStoichiometries));
            }
            testProt = testProt.Where(p => p.ProteinPairwiseComparisons != null).ToList();
            Extensions.CalcCorrectedPValue(testProt, groupPepsForPValCalc, alpha);
            List<ProteinGroup> ProteinsToUse = testProt.Where(p => p.useProt).ToList();
            //store protein output in XML file
            string outFile = @"D:\PTMStoichiometry\TestData\EcoliSpikeIn\2021-06-03-11-06-06\Task3-SearchTask\20210604a_AllQuantifiedPeptidesEvenGroupsPeptidePeptideUngrouped.xml"; //replace with desired output file
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
