using System;
using System.Collections.Generic;
using System.Text;

namespace PTMStoichiometry20210414a
{
    //this class calculates the stoichmetry of a peptide in a group: pep/baseline
    public class Stoichiometry
    {
        public double StoichiometryVal { get;  } //intensity of pep/baseline or pep1/pep2 ;
        public Boolean usefulStoichiometry { get; } //true if pep Intensity > 0, false otherwise
        public Stoichiometry(Intensity Pep, double baseline)
        { 
            this.StoichiometryVal = Pep.IntensityVal / baseline;
            this.usefulStoichiometry = (Pep.IntensityVal > 0); 
        }

        //overload for peptide:peptide case
        public Stoichiometry(Intensity Pep1, Intensity Pep2)
        {
            this.StoichiometryVal = Pep1.IntensityVal / Pep2.IntensityVal;
            this.usefulStoichiometry = (Pep1.IntensityVal > 0 && Pep2.IntensityVal > 0);
        }

    }
}
