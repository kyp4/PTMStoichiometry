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
            int reqNumModPeptides = 3;
            int reqNumOfPepeptides = reqNumUnmodPeptides + reqNumModPeptides;
            bool useBaselinePeptides = true;
            int reqNumBaselinePeptides = reqNumUnmodPeptides;
            int reqNumBaselineMeasurements = 5; //allow one missing value
            double correlationCutOff = 0.5;
            bool compareUnmod = false;
            int minNumStoichiometries = 3;
            int reqNumPepMeasurements = 3;
            bool groupPepsForPValCalc = true;
            double alpha = 0.05;
            string groupToCompare = null;
            string dataType = "FlashLFQ";

            string filepathpeptides = @"D:\PTMStoichiometry\TestData\MSV000086126\2021-06-03-16-11-24\Task3-SearchTask\AllQuantifiedPeptides.tsv";
            string filepathgroups = @"D:\PTMStoichiometry\TestData\MSV000086126\2021-06-03-16-11-24\Task3-SearchTask\MSV000086126Groups.txt";
            string directory = @"D:\PTMStoichiometry\TestData\MSV000086126\2021-06-03-16-11-24\Task3-SearchTask";
            string stoichiometryfileout = "20210607g_MSV000086126";
            string paramsfile = stoichiometryfileout + "_params";
            
            Dictionary<string, string> groups = PeptideReader.GetGroups(filepathgroups);
            List<string> groupsList = PeptideReader.GetGroupList(filepathgroups);
            List<Peptide> testPeptide = PeptideReader.ReadTsv(filepathpeptides, filepathgroups, reqNumPepMeasurements, dataType);
            testPeptide = Extensions.IncludeSharedPeptides(testPeptide, false); 


            //group peptides by protein
            var proteins = testPeptide.Select(p => p.ProteinGroup).Distinct().ToArray();
            List<ProteinGroup> testProt = new List<ProteinGroup>();

            if (groupToCompare != null)
            {
                for (int i = 0; i < proteins.Length; i++)
                {
                    testProt.Add(new ProteinGroup(proteins[i], testPeptide, groupsList, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides,
                    useBaselinePeptides, reqNumBaselinePeptides, reqNumBaselineMeasurements, correlationCutOff, compareUnmod, minNumStoichiometries, groupToCompare));
                }
            }
            else
            {
                for (int i = 0; i < proteins.Length; i++)
                {
                    testProt.Add(new ProteinGroup(proteins[i], testPeptide, groupsList, reqNumUnmodPeptides, reqNumModPeptides, reqNumOfPepeptides,
                    useBaselinePeptides, reqNumBaselinePeptides, reqNumBaselineMeasurements, correlationCutOff, compareUnmod, minNumStoichiometries));
                }
            }

            testProt = testProt.Where(p => p.ProteinPairwiseComparisons != null).ToList(); //hmmm
            Extensions.CalcCorrectedPValue(testProt, groupPepsForPValCalc, alpha);
            List<ProteinGroup> ProteinsToUse = testProt.Where(p => p.useProt).ToList();
            
            WriteFile.ParamsWriter(paramsfile, filepathpeptides, filepathgroups, directory, stoichiometryfileout, reqNumUnmodPeptides, reqNumModPeptides,
                reqNumOfPepeptides, useBaselinePeptides, reqNumBaselinePeptides, reqNumBaselineMeasurements, correlationCutOff, compareUnmod,
                minNumStoichiometries, reqNumPepMeasurements, groupPepsForPValCalc, alpha, groupToCompare);
            WriteFile.StoichiometryDataWriter(ProteinsToUse, useBaselinePeptides, minNumStoichiometries, directory, stoichiometryfileout);

        }
    }
}
