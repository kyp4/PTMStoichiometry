using NUnit.Framework;
using System;
using System.Collections.Generic;
using PTMStoichiometry;
using System.Linq;
using System.IO;

namespace Test
{
    [TestFixture]
    class Class1
    {

        
        private static readonly object[] _calcTest =
        {
            new object[] {
                "KTEM[Common Variable:Oxidation on M]VSSVPAE[Metal:FeIII on E]NKSVLNEHQETSK",
                 new List<string>() { "Common Variable:Oxidation on M", "Metal:FeIII on E" },
                 "FlashLFQ",
                 new List<string>() { "KTEMCommon Variable:Oxidation on MVSSVPAENKSVLNEHQETSK", "KTEMVSSVPAEMetal:FeIII on ENKSVLNEHQETSK" }
            },
            new object[] {
                "[Common Biological:Acetylation on X]RAEEPC[Common Fixed:Carbamidomethyl on C]APGAPSALGAQR",
                 new List<string>() { "Common Biological:Acetylation on X", "Common Fixed:Carbamidomethyl on C" },
                 "FlashLFQ",
                 new List<string>() { "RAEEPCCommon Fixed:Carbamidomethyl on CAPGAPSALGAQR" }
            },
        };

        [Test]
        //KTEM[Common Variable:Oxidation on M]VSSVPAE[Metal:Fe[III] on E]NKSVLNEHQETSK
        //[Common Biological:Acetylation on X]RAEEPC[Common Fixed:Carbamidomethyl on C]APGAPSALGAQR
        [TestCaseSource("_calcTest")]
        public void Peptide_GetLocalizedModifications_Pass(string seq, List<string> mods, string dataType, List<string> localized)
        {
            List<string> testGetLocalizedModifications = Peptide.GetLocalizedModifications(seq, mods, dataType);
    
            Assert.AreEqual(localized, testGetLocalizedModifications);
        }

    }
}
