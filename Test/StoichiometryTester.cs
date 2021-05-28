using NUnit.Framework;
using System;
using System.Collections.Generic;
using PTMStoichiometry20210414a;
using System.Linq;

namespace PTMStoichiometryTester20200415a
{
    [TestFixture]
    class StoichiometryTester
    {

        [Test]
        [TestCase("file.txt", "Control", 284.6, DetectionMS.MSMS, 45737.7547, true)]
        [TestCase("file.txt", "Control", 2654754284.5686586, DetectionMS.MSMS, 45737.7547, true)]
        [TestCase("file3.txt", "group 45", 3579.98456, DetectionMS.MS, 27346.3456, true)]
        [TestCase("file3.txt", "group 45", 9843.43, DetectionMS.MS, 273443536.3456, true)]
        [TestCase("file3.raw", "group-42%AA", 0, DetectionMS.MSMSIdentifiedButNotQuantified, 3857.48, false)]
        [TestCase("file 4", "98%HA", 0, DetectionMS.NotDetected, 476.64, false)]
        [TestCase("file_46_4.raw", "Disease 01", 0, DetectionMS.MSMSAmbiguousPeakfinding, 1461636.42, false)]
        public void Stoichiometry_BaselineCase_Pass(string file, string group, double intensity, DetectionMS detection, double baseline, Boolean useful)
        {
            Intensity Testintensity = new Intensity(file, group, intensity, detection);
            Stoichiometry TestStoichiometry = new Stoichiometry(Testintensity, baseline);
            Assert.AreEqual(Testintensity.IntensityVal / baseline, TestStoichiometry.StoichiometryVal);
            Assert.AreEqual(useful, TestStoichiometry.usefulStoichiometry);
        }

        [Test]
        [TestCase("file.txt", "Control", 284.6, 82746.347, DetectionMS.MSMS, DetectionMS.MS, true)]
        [TestCase("file.txt", "Control", 2654754284.5686586, 45737.7547, DetectionMS.MSMS, DetectionMS.MSMS, true)]
        [TestCase("file3.txt", "group 45", 3579.98456, 27346.3456, DetectionMS.MS, DetectionMS.MS, true)]
        [TestCase("file3.txt", "group 45", 9843.43, 273443536.3456, DetectionMS.MS,  true)]
        [TestCase("file3.raw", "group-42%AA", 0, 3857.48, DetectionMS.MSMSIdentifiedButNotQuantified, DetectionMS.MS, false)]
        [TestCase("file 4", "98%HA", 476.64, 0, DetectionMS.MSMS, DetectionMS.NotDetected, false)]
        [TestCase("file_46_4.raw", "Disease 01", 0, 0, DetectionMS.MSMSAmbiguousPeakfinding, DetectionMS.MSMSIdentifiedButNotQuantified, false)]
        public void Stoichiometry_PeptidePeptideCase_Pass(string file, string group, double intensity1, double intensity2, DetectionMS detection1, DetectionMS detection2, Boolean useful)
        {
            Intensity Testintensity1 = new Intensity(file, group, intensity1, detection1);
            Intensity Testintensity2 = new Intensity(file, group, intensity2, detection2);
            Stoichiometry TestStoichiometry = new Stoichiometry(Testintensity1, Testintensity2);
            Assert.AreEqual(Testintensity1.IntensityVal / Testintensity2.IntensityVal, TestStoichiometry.StoichiometryVal);
            Assert.AreEqual(useful, TestStoichiometry.usefulStoichiometry);
        }

    }
}
