using System;
using System.Collections.Generic;
using System.Text;

namespace PTMStoichiometry20210414a
{
    //this class calculates the stoichmetry of a peptide in a group: pep/baseline
    public class Stoichiometry
    {
        public double StoichiometryVal { get;  } //intensity of pep/baseline
        public Boolean usefulStoichiometry { get; } //true if pep Intensity is MS or MSMS (so has values), false otherwise
        public Stoichiometry(Intensity Pep, double baseline)
        { 
            this.StoichiometryVal = Pep.IntensityVal / baseline;
            this.usefulStoichiometry = ((Pep.Detection == DetectionMS.MS || Pep.Detection == DetectionMS.MSMS)); //check that there is intensity to compare
        }

    }
}
