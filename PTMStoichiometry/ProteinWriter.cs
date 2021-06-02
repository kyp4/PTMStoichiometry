using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.Linq;

namespace PTMStoichiometry
{
    //class to write protein output to XML file TODO: need dif writer for peptide and baseline cases
    public class ProteinWriter
    {
        //code to add protein node to Xml file: adapted from https://docs.microsoft.com/en-us/dotnet/api/system.xml.xmldocument?view=net-5.0
        public static XmlElement AddProtein(ProteinGroup prot, XmlDocument doc)
        { 
            //create a node for each protein group
            XmlElement protNode = doc.CreateElement("ProteinGroup");

            //protein group attributes: accesion, number of peptides in protein, Benjamini-Hochberg corrected p-value for protein
            XmlAttribute attribute = doc.CreateAttribute("ProteinAcession");
            attribute.Value = prot.ProteinName;
            protNode.Attributes.Append(attribute);

            attribute = doc.CreateAttribute("NumberOfPepetides");
            attribute.Value = prot.NumPeptidesInProtein.ToString();
            protNode.Attributes.Append(attribute);

            //attribute = doc.CreateAttribute("BaselinePeptides");
            //attribute.Value = String.Join(";", prot.BaselinePeptides.Select(p => p.Sequence));
            //protNode.Attributes.Append(attribute);

            //Create child node for each pairwise comparison on peptides 
            foreach (PairwiseCompairison comp in prot.ProteinPairwiseComparisons)
            {
                if (comp.PeptideStoichiometriesGroupOne.Count() > 3 && comp.PeptideStoichiometriesGroupTwo.Count() > 3) //only include peptides with enough data to calc stats
                {
                    XmlElement peptideComp = doc.CreateElement("PairwiseCompairison");
                    //add pairwise compairison attributes: Peptides compaired, groups compared
                    XmlAttribute attributePep = doc.CreateAttribute("Peptide");
                    attributePep.Value = comp.PeptideOne.Sequence;
                    peptideComp.Attributes.Append(attributePep);

                    attributePep = doc.CreateAttribute("Groups");
                    attributePep.Value = comp.GroupOne + " | " + comp.GroupTwo;
                    peptideComp.Attributes.Append(attributePep);

                    //add pairwise comparison info: statistics, stoichiometries as children nodes
                    XmlElement MWStat = doc.CreateElement("Mann-Whitney-TestStat");
                    MWStat.InnerText = comp.MWStat.ToString();
                    peptideComp.AppendChild(MWStat);

                    XmlElement MWpval = doc.CreateElement("Mann-Whitney-p-value");
                    MWpval.InnerText = comp.MWPVal.ToString();
                    peptideComp.AppendChild(MWpval);

                    XmlElement bhCorrectedpvalue = doc.CreateElement("Benjamini-Hochberg-corrected-p-value");
                    bhCorrectedpvalue.InnerText = comp.CorrectedpVal.ToString();
                    peptideComp.AppendChild(bhCorrectedpvalue);

                    XmlElement StoicMedG1 = doc.CreateElement("MedianStoichiometryGroupOne");
                    StoicMedG1.InnerText = comp.PeptideStoichiometriesGroupOneMedian.ToString();
                    peptideComp.AppendChild(StoicMedG1);

                    XmlElement StoicMinG1 = doc.CreateElement("MinStoichiometryGroupOne");
                    StoicMinG1.InnerText = comp.PeptideStoichiometriesGroupOneMin.ToString();
                    peptideComp.AppendChild(StoicMinG1);

                    XmlElement StoicMaxG1 = doc.CreateElement("MaxStoichiometryGroupOne");
                    StoicMaxG1.InnerText = comp.PeptideStoichiometriesGroupOneMax.ToString();
                    peptideComp.AppendChild(StoicMaxG1);

                    XmlElement StoicMedG2 = doc.CreateElement("MedianStoichiometryGroupTwo");
                    StoicMedG2.InnerText = comp.PeptideStoichiometriesGroupTwoMedian.ToString();
                    peptideComp.AppendChild(StoicMedG2);

                    XmlElement StoicMinG2 = doc.CreateElement("MinStoichiometryGroupTwo");
                    StoicMinG2.InnerText = comp.PeptideStoichiometriesGroupTwoMin.ToString();
                    peptideComp.AppendChild(StoicMinG2);

                    XmlElement StoicMaxG2 = doc.CreateElement("MaxStoichiometryGroupTwo");
                    StoicMaxG2.InnerText = comp.PeptideStoichiometriesGroupTwoMax.ToString();
                    peptideComp.AppendChild(StoicMaxG2);

                    /*
                    foreach (Stoichiometry s1 in comp.PeptideStoichiometriesGroupOne)
                    {
                        XmlElement StoichG1 = doc.CreateElement("StoichiometryGroupOne");
                        StoichG1.InnerText = s1.StoichiometryVal.ToString();
                        peptideComp.AppendChild(StoichG1);
                    }

                    foreach (Stoichiometry s2 in comp.PeptideStoichiometriesGroupTwo)
                    {
                        XmlElement StoichG2 = doc.CreateElement("StoichiometryGroupTwo");
                        StoichG2.InnerText = s2.StoichiometryVal.ToString();
                        peptideComp.AppendChild(StoichG2);
                    }
                    */
                    protNode.AppendChild(peptideComp);
                }
            }

            return protNode;
        }

    }
}
