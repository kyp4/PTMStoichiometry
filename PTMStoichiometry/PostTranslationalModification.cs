using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace PTMStoichiometry
{
    //class to store post translational modification data
    public class PostTranslationalModification
    {
        //whether the PTM is localized 
        public bool Localized { get; }
        //the modificationon it's own
        public string Modification { get; }
        //the modification within the peptide sequence with all other ptms removed
        public string ModificationInPeptideSequence { get; }
        //peptide seq 
        public string Peptide { get; }

        /// <summary>
        /// Object to hold post translational modification data
        /// </summary>
        /// <param name="modification">the modificationon it's own</param>
        /// <param name="modificationInPeptideSequence">the modification within the peptide sequence with 
        /// all other ptms removed</param>
        public PostTranslationalModification(string modification, string modificationInPeptideSequence, string seq)
        {
            this.Modification = modification;
            this.Localized = !this.Modification.Contains("on X");
            this.ModificationInPeptideSequence = modificationInPeptideSequence;
            this.Peptide = seq;
        }

    }
}
