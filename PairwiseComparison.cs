using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Meta.Numerics.Statistics;

namespace PTMStoichiometry20210414a
{
    //class for performing a pairwise compairson of peptide data between two groups
    public class PairwiseComparison
    {
        //Peptides being compared
        public Peptide PeptideOne { get; }
        public Peptide PeptideTwo { get; }

        //groups being compared
        public string GroupOne { get; }
        public string GroupTwo { get; }
        //peptide Stoichiometries
        public List<Stoichiometry> PeptideStoichiometriesGroupOne { get; }
        public List<Stoichiometry> PeptideStoichiometriesGroupTwo { get; }

        //peptide Stoichiometries - median and min/max for each group
        public double PeptideStoichiometriesGroupOneMedian { get; }
        public double PeptideStoichiometriesGroupTwoMedian { get; }
        public double PeptideStoichiometriesGroupOneMin { get; }
        public double PeptideStoichiometriesGroupTwoMin { get; }
        public double PeptideStoichiometriesGroupOneMax { get; }
        public double PeptideStoichiometriesGroupTwoMax { get; }
        //Mann-Whitney test stat comparing stoichiometries between the two groups 
        public double MWStat { get; }
        //p-value
        public double MWPVal { get; }

        public PairwiseComparison(Peptide Pep1, Peptide Pep2, string G1, string G2)
        {
            this.PeptideOne = Pep1;
            this.PeptideTwo = Pep2;
            this.GroupOne = G1;
            this.GroupTwo = G2;
            this.PeptideStoichiometriesGroupOne = calcStoichiometry(this.GroupOne);
            this.PeptideStoichiometriesGroupTwo = calcStoichiometry(this.GroupTwo);
            //private double[] mwStats = calcMWStats();
            if (this.PeptideStoichiometriesGroupOne.Count > 3 && this.PeptideStoichiometriesGroupTwo.Count > 3) //require at least three stoichiometries in both dist to run test
            {
                this.MWStat = calcMWStats()[0];
                this.MWPVal = calcMWStats()[1];
                this.PeptideStoichiometriesGroupOneMedian = calcMWStats()[2];
                this.PeptideStoichiometriesGroupTwoMedian = calcMWStats()[3];
                this.PeptideStoichiometriesGroupOneMin = calcMWStats()[4];
                this.PeptideStoichiometriesGroupTwoMin = calcMWStats()[5];
                this.PeptideStoichiometriesGroupOneMax = calcMWStats()[6];
                this.PeptideStoichiometriesGroupTwoMax = calcMWStats()[7];
            }
        }

        //calculate stoichiometries for all intensities
        private List<Stoichiometry> calcStoichiometry(string group)
        {
            List<Stoichiometry> stoich = new List<Stoichiometry>();
            List<Intensity> Pep1Intensity = this.PeptideOne.Intensities.Where(p => p.GroupID == group).ToList(); //intensities pep1 for group of interest
            List<Intensity> Pep2Intensity = this.PeptideTwo.Intensities.Where(p => p.GroupID == group).ToList(); //intensities pep2 for group of interest

            foreach (Intensity i1 in Pep1Intensity)
            {
                foreach (Intensity i2 in Pep2Intensity)
                {
                    Stoichiometry calcStoich = new Stoichiometry(i1, i2);
                    if (calcStoich.usefulStoichiometry)
                    {
                        stoich.Add(new Stoichiometry(i1, i2));
                    }
                }
            }
            return stoich;
        }

        //function to calc Mann-Whitney statistics
        private double[] calcMWStats()
        {
            double[] stats = new double[8]; //to store final stats
            List<double> stoichsG1 = new List<double>(); //double[(this.PeptideStoichiometriesGroupOne.Count()]; //list to store distribution of group one stoichiometries
            List<double> stoichsG2 = new List<double>(); //double[(this.PeptideStoichiometriesGroupTwo.Count()]; //list to store distribution of group two stoichiometries
            foreach (Stoichiometry S1 in (this.PeptideStoichiometriesGroupOne)) {
                stoichsG1.Add(S1.StoichiometryVal);
            }
            foreach (Stoichiometry S2 in (this.PeptideStoichiometriesGroupTwo)) {
                stoichsG2.Add(S2.StoichiometryVal);
            }
            //calc stats
            TestResult mw = Univariate.MannWhitneyTest(stoichsG1, stoichsG2); 
            stats[0] = mw.Statistic.Value;
            stats[1] = mw.Probability;

            //find medians, mins maxs
            stats[2] = stoichsG1.Median(); 
            stats[4] = stoichsG1.Min();
            stats[6] = stoichsG1.Max();

            stats[3] = stoichsG2.Median();
            stats[5] = stoichsG2.Min();
            stats[7] = stoichsG2.Max();

            return stats;
        }
    }
}
