using System;
using System.Collections.Generic;
using System.Text;

namespace PTMStoichiometry20210414a
{
    //this class calculates the stoichmetry between two peptides in a group: pep1/pep2
    public class Stoichiometry
    {
        public double StoichiometryVal { get;  } //intensity of pep1/pep2
        public Boolean usefulStoichiometry { get; } //true if both Intensities are MS or MSMS (so have values), false otherwise
        public Stoichiometry(Intensity Pep1, Intensity Pep2)
        { 
            this.StoichiometryVal = Pep1.IntensityVal / Pep2.IntensityVal;
            this.usefulStoichiometry = ((Pep1.Detection == DetectionMS.MS || Pep1.Detection == DetectionMS.MSMS) &&
                                            (Pep2.Detection == DetectionMS.MS || Pep2.Detection == DetectionMS.MSMS) &&
                                            (Pep1 != Pep2)); //check that there is intensity to compare and the Peptide is not being compared to self
        }

    }
}
