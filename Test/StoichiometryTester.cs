using NUnit.Framework;
using System;
using System.Collections.Generic;
using PTMStoichiometry;
using System.Linq;

namespace Test
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
        [TestCase("file.txt",       "Control", 284.6, 82746.347, DetectionMS.MSMS, DetectionMS.MS, true, 0.003439427)]
        [TestCase("file.txt",       "Control", 2654754284.5686586, 45737.7547, DetectionMS.MSMS, DetectionMS.MSMS, true, 58042.952)]
        [TestCase("file3.txt",      "group 45", 3579.98456, 27346.3456, DetectionMS.MS, DetectionMS.MS, true, 0.1309127)]
        [TestCase("file3.txt",      "group 45", 9843.43, 273443536.3456, DetectionMS.MS, DetectionMS.MS, true, 0.00003599804)]
        [TestCase("file3.raw",      "group-42%AA", 0, 3857.48, DetectionMS.MSMSIdentifiedButNotQuantified, DetectionMS.MS, false, 0)]
        [TestCase("file 4",         "98%HA",     476.64, 0, DetectionMS.MSMS, DetectionMS.NotDetected, false, double.PositiveInfinity)]
        [TestCase("file_46_4.raw",  "Disease 01", 0, 0, DetectionMS.MSMSAmbiguousPeakfinding, DetectionMS.MSMSIdentifiedButNotQuantified, false, double.NaN)]
        public void Stoichiometry_PeptidePeptideCase_Pass(string file, string group, double intensity1, double intensity2, DetectionMS detection1, 
            DetectionMS detection2, bool useful, double stoich)
        {
            Intensity Testintensity1 = new Intensity(file, group, intensity1, detection1);
            Intensity Testintensity2 = new Intensity(file, group, intensity2, detection2);
            Stoichiometry TestStoichiometry = new Stoichiometry(Testintensity1, Testintensity2);
            Assert.That(TestStoichiometry.StoichiometryVal, Is.EqualTo(stoich).Within(0.001));
            Assert.AreEqual(useful, TestStoichiometry.usefulStoichiometry);
        }
    }
}
