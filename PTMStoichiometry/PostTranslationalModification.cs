using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace PTMStoichiometry
{
    public class PostTranslationalModification
    {
        public bool Localized { get; }
        public string Modification { get; }
        public string ModificationInPeptideSequence { get; }

        public PostTranslationalModification(string modification, string modificationInPeptideSequence)
        {
            this.Modification = modification;
            this.Localized = !this.Modification.Contains("on X");
            this.ModificationInPeptideSequence = modificationInPeptideSequence;
        }

    }
}
