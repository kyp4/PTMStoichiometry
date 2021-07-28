using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;

namespace PTMStoichiometry
{
    //this class calculates the stoichmetry of a peptide in a group
    public class Stoichiometry
    {
        //intensity of pep(s)/baseline or pep1/pep2
        public double StoichiometryVal { get; }
        //true if pep Intensity > 0, false otherwise
        public bool usefulStoichiometry { get; } 
       
        /// <summary>
        /// object to store stoichiometry dat: peptides (same ptm)/baseline
        /// </summary>
        /// <param name="Pep">list of intensities from a peptide</param>
        /// <param name="baseline">baseline value - the denominator</param>
        public Stoichiometry(List<Intensity> Pep, double baseline)
        { 
            this.StoichiometryVal = Pep.Select(p => p.IntensityVal).Sum() / baseline;
            this.usefulStoichiometry = (Pep.Select(p => p.IntensityVal).Sum() > 0); 
        }

        /// <summary>
        /// object to store stoichiometry dat: peptide/baseline
        /// </summary>
        /// <param name="Pep">single intensity from a peptide</param>
        /// <param name="baseline">baseline value - the denominator</param>
        public Stoichiometry(Intensity Pep, double baseline)
        {
            this.StoichiometryVal = Pep.IntensityVal / baseline;
            this.usefulStoichiometry = (Pep.IntensityVal > 0);
        }

        /// <summary>
        /// object to store stoichiometry data: peptides (same ptm)/peptide
        /// </summary>
        /// <param name="Pep1">list of intensities from a peptide</param>
        /// <param name="Pep2">single intensity from a peptide - the denominator</param>
        public Stoichiometry(List<Intensity> Pep1, Intensity Pep2)
        {
            this.StoichiometryVal = Pep1.Select(p => p.IntensityVal).Sum() / Pep2.IntensityVal;
            this.usefulStoichiometry = (Pep1.Select(p => p.IntensityVal).Sum() > 0 && Pep2.IntensityVal > 0);
        }


        /// <summary>
        /// object to store stoichiometry data: peptide/peptide
        /// </summary>
        /// <param name="Pep1">single intensity from a peptide</param>
        /// <param name="Pep2">single intensity from a peptide - the denominator</param>
        public Stoichiometry(Intensity Pep1, Intensity Pep2)
        {
            this.StoichiometryVal = Pep1.IntensityVal / Pep2.IntensityVal;
            this.usefulStoichiometry = (Pep1.IntensityVal > 0 && Pep2.IntensityVal > 0);
        }

    }
}
