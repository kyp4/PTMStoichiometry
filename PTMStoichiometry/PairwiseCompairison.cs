using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Statistics;
using Meta.Numerics.Statistics;

namespace PTMStoichiometry
{
    //class for performing a pairwise compairson of peptide data between two groups
    public class PairwiseCompairison
    {
        //Peptides being compared (in baseline case, there is only one peptide)
        public Peptide PeptideOne { get; }
        public Peptide PeptideTwo { get; }

        public List<Peptide> PeptidesWithPTM { get; }
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
        public double CorrectedpVal { get; set; }
        public string PTM { get; }

        //minNumStoichiometries - min num of stoichiometries req in both groups before run test (default=3)
        public PairwiseCompairison(Peptide Pep, List<Intensity> BaselinePepsIntensity, string G1, string G2, int minNumStoichiometries)
        {
            this.PeptideOne = Pep;
            this.GroupOne = G1;
            this.GroupTwo = G2;
            this.PeptideStoichiometriesGroupOne = calcStoichiometry(this.PeptideOne, this.GroupOne, BaselinePepsIntensity);
            this.PeptideStoichiometriesGroupTwo = calcStoichiometry(this.PeptideOne, this.GroupTwo, BaselinePepsIntensity);
            
            if (this.PeptideStoichiometriesGroupOne.Count() >= minNumStoichiometries && this.PeptideStoichiometriesGroupTwo.Count() >= minNumStoichiometries) 
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

        public PairwiseCompairison(List<Peptide> Peps, List<Intensity> BaselinePepsIntensity, string G1, string G2, int minNumStoichiometries, string ptm)
        {
            this.PeptidesWithPTM = Peps;
            this.PTM = ptm;
            this.GroupOne = G1;
            this.GroupTwo = G2;
            this.PeptideStoichiometriesGroupOne = calcStoichiometry(this.PeptidesWithPTM, this.GroupOne, BaselinePepsIntensity);
            this.PeptideStoichiometriesGroupTwo = calcStoichiometry(this.PeptidesWithPTM, this.GroupTwo, BaselinePepsIntensity);

            if (this.PeptideStoichiometriesGroupOne.Count() >= minNumStoichiometries && this.PeptideStoichiometriesGroupTwo.Count() >= minNumStoichiometries)
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
        //overload: for peptide to peptide case
        public PairwiseCompairison(Peptide Pep1, Peptide Pep2, string G1, string G2, int minNumStoichiometries)
        {
            this.PeptideOne = Pep1;
            this.PeptideTwo = Pep2;
            this.GroupOne = G1;
            this.GroupTwo = G2;
            this.PeptideStoichiometriesGroupOne = calcStoichiometry(this.PeptideOne, this.PeptideTwo, this.GroupOne);
            this.PeptideStoichiometriesGroupTwo = calcStoichiometry(this.PeptideOne, this.PeptideTwo, this.GroupTwo);

            if (this.PeptideStoichiometriesGroupOne.Select( p => p.usefulStoichiometry).Count() >= minNumStoichiometries 
                && this.PeptideStoichiometriesGroupTwo.Select(p => p.usefulStoichiometry).Count() >= minNumStoichiometries)
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

        public void setCorrectedpVal(double correctPVal)
        {
            this.CorrectedpVal = correctPVal;
        }

        private List<Stoichiometry> calcStoichiometry(Peptide pep, string group, List<Intensity> baselineIntensity)
        {
            List<Stoichiometry> stoich = new List<Stoichiometry>();
            List<Intensity> PepIntensity = pep.Intensities.Where(p => p.GroupID == group).ToList(); //intensities pep1 for group of interest
            List<Intensity> baselineGroupIntensity = baselineIntensity.Where(p => p.GroupID == group).ToList();

            double baseline = baselineGroupIntensity.Select(p => p.IntensityVal).Median();
            for (int i = 0; i < PepIntensity.Count(); i++)
            {                
                stoich.Add(new Stoichiometry(PepIntensity[i], baseline));
            }
            return stoich;
        }

        //calculate stoichiometries for all intensities
        private List<Stoichiometry> calcStoichiometry(List<Peptide> peps, string group, List<Intensity> baselineIntensity)
        {
            List<Stoichiometry> stoich = new List<Stoichiometry>();
            //List<Intensity> PepIntensity = peps.Select(p => p.Intensities.Where(p => p.GroupID == group)).ToList(); //intensities pep1 for group of interest
            List<Intensity> baselineGroupIntensity = baselineIntensity.Where(p => p.GroupID == group).ToList();

            List<string> baselineFileNames = baselineGroupIntensity.Select(p => p.FileName).Distinct().ToList();

            double baseline = baselineGroupIntensity.Select(p => p.IntensityVal).Median();
            for (int i = 0; i < baselineFileNames.Count(); i++)
            {
                List<Intensity> PepIntensity = new List<Intensity>();
                foreach (Peptide pep in peps)
                {
                    //List<Intensity> pepInt = pep.Intensities.Where(p => p.GroupID == group).ToList();
                    PepIntensity.Add(pep.Intensities.Where(p => p.FileName.Equals(baselineFileNames[i])).ToList()[0]);
                }
                stoich.Add(new Stoichiometry(PepIntensity, baseline));               
            }
            return stoich;
        }

        //overload: peptide to peptide case
        private List<Stoichiometry> calcStoichiometry(Peptide pep1, Peptide pep2, string group)
        {
            List<Stoichiometry> stoich = new List<Stoichiometry>();
            List<Intensity> Pep1Intensity = pep1.Intensities.Where(p => p.GroupID == group).ToList();
            List<Intensity> Pep2Intensity = pep2.Intensities.Where(p => p.GroupID == group).ToList();

            foreach (Intensity i1 in Pep1Intensity)
            {
                Intensity i2 = Pep2Intensity.Where(p => p.FileName == i1.FileName).ToList()[0];
                Stoichiometry calcStoich = new Stoichiometry(i1, i2);
                if (i1.IntensityVal > 0 || i2.IntensityVal > 0)
                {
                    stoich.Add(calcStoich);
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
