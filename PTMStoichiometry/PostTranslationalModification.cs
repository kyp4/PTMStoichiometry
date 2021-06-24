using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace PTMStoichiometry
{
    class PostTranslationalModification
    {
        bool Localized { get; }
        string Modification { get; }
        string ModificationInPeptideSequence { get; }

        public PostTranslationalModification(string modification, string modificationInPeptideSequence)
        {
            this.Modification = modification;
            this.Localized = !this.Modification.Contains("in X");
            this.ModificationInPeptideSequence = modificationInPeptideSequence;


        }

    }
}
