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
        /// <summary>
        /// Test to check that stoichiometry is being calculated correctly in the for a peptide with a baseline
        /// </summary>
        /// <param name="file">the filethe intensity data came from</param>
        /// <param name="group">the group the intensity is associated with</param>
        /// <param name="intensity">the intensity value</param>
        /// <param name="baseline">the denominator of the stoichiometry equation</param>
        /// <param name="useful">correct bool to compare to usefulStoichiometry</param>
        [Test]
        [TestCase("file.txt", "Control", 284.6, 45737.7547, true)]
        [TestCase("file.txt", "Control", 2654754284.5686586, 45737.7547, true)]
        [TestCase("file3.txt", "group 45", 3579.98456, 27346.3456, true)]
        [TestCase("file3.txt", "group 45", 9843.43, 273443536.3456, true)]
        [TestCase("file3.raw", "group-42%AA", 0, 3857.48, false)]
        [TestCase("file 4", "98%HA", 0, 476.64, false)]
        [TestCase("file_46_4.raw", "Disease 01", 0,  1461636.42, false)]
        public void Stoichiometry_PeptideBaselineCase_Pass(string file, string group, double intensity, 
            double baseline, bool useful)
        {
            Intensity Testintensity = new Intensity(file, group, intensity);
            Stoichiometry TestStoichiometry = new Stoichiometry(Testintensity, baseline);
            Assert.AreEqual(Testintensity.IntensityVal / baseline, TestStoichiometry.StoichiometryVal);
            Assert.AreEqual(useful, TestStoichiometry.usefulStoichiometry);
        }

        /// <summary>
        /// Test to check that stoichiometry is being calculated correctly in the for a peptide in the peptide:peptide case
        /// </summary>
        /// <param name="file">the filethe intensity data came from</param>
        /// <param name="group">the group the intensity is associated with</param>
        /// <param name="intensity1">the intensity value of the numerator</param>
        /// <param name="intensity2">the intensity value of the denominator</param>
        /// <param name="useful">correct bool to compare to usefulStoichiometry</param>
        [Test]
        [TestCase("file.txt",       "Control",      284.6,              82746.347,       true)]
        [TestCase("file.txt",       "Control",      2654754284.5686586, 45737.7547,      true)]
        [TestCase("file3.txt",      "group 45",     3579.98456,         27346.3456,      true)]
        [TestCase("file3.txt",      "group 45",     9843.43,            273443536.3456,  true)]
        [TestCase("file3.raw",      "group-42%AA",  0,                  3857.48,         false)]
        [TestCase("file 4",         "98%HA",        476.64,             0,               false)]
        [TestCase("file_46_4.raw",  "Disease 01",   0,                  0,               false)]
        public void Stoichiometry_PeptidePeptidePeptideCase_Pass(string file, string group, double intensity1,
            double intensity2, bool useful)
        {
            Intensity Testintensity1 = new Intensity(file, group, intensity1);
            Intensity Testintensity2 = new Intensity(file, group, intensity2);
            Stoichiometry TestStoichiometry = new Stoichiometry(Testintensity1, Testintensity2);
            Assert.That(TestStoichiometry.StoichiometryVal, Is.EqualTo(intensity1 / intensity2).Within(0.001));
            Assert.AreEqual(useful, TestStoichiometry.usefulStoichiometry);
        }


        private static readonly object[] _testPTMBaselineCase =
        {
            new object[] {
                "file.txt",
                "Control",
                new List<double> { 284.6 },
                45737.7547,
                true
            },
            new object[] {
                "file.txt",
                "Control",
                new List<double> { 2654754284.5686586, 3579.98456 },
                45737.7547,
                true
            },
            new object[] {
                "file.txt",
                "Control",
                new List<double> { 0, 3579.98456 },
                45737.7547,
                true
            },
            new object[] {
                "file.txt",
                "Control",
                new List<double> { 0 },
                45737.7547,
                false
            },
            new object[] {
                "file.txt",
                "Control",
                new List<double> { 2654754284.5686586, 3579.98456, 0, 0, 89353.9834, 0, 239847 },
                45737.7547,
                true
            },
            new object[] {
                "file 4", 
                "98%HA",
                new List<double> { 0, 0 },
                45737.7547,
                false
            }
        };
        /// <summary>
        /// Test to check that stoichiometry is being calculated correctly in the for peptides with the same ptm with a baseline
        /// </summary>
        /// <param name="file">the filethe intensity data came from</param>
        /// <param name="group">the group the intensity is associated with</param>
        /// <param name="intensities">the intensities values</param>
        /// <param name="baseline">the denominator of the stoichiometry equation</param>
        /// <param name="useful">correct bool to compare to usefulStoichiometry</param>
        [Test]
        [TestCaseSource("_testPTMBaselineCase")]
        public void Stoichiometry_PTMBaselineCase_Pass(string file, string group, List<double> intensities,
            double baseline, bool useful)
        {
            List<Intensity> Testintensities = new List<Intensity>();
            foreach (double d in intensities)
            {
                Testintensities.Add(new Intensity(file, group, d));
            }

            Stoichiometry TestStoichiometry = new Stoichiometry(Testintensities, baseline);
            Assert.AreEqual(intensities.Sum() / baseline, TestStoichiometry.StoichiometryVal);
            Assert.AreEqual(useful, TestStoichiometry.usefulStoichiometry);
        }

        private static readonly object[] _testPTMPeptidePeptideCase =
        {
            new object[] {
                "file.txt",
                "Control",
                new List<double> { 284.6 },
                45737.7547,
                true
            },
            new object[] {
                "file.txt",
                "Control",
                new List<double> { 2654754284.5686586, 3579.98456 },
                45737.7547,
                true
            },
            new object[] {
                "file.txt",
                "Control",
                new List<double> { 0, 3579.98456 },
                45737.7547,
                true
            },
            new object[] {
                "file.txt",
                "Control",
                new List<double> { 0 },
                45737.7547,
                false
            },
            new object[] {
                "file.txt",
                "Control",
                new List<double> { 2654754284.5686586, 3579.98456, 0, 0, 89353.9834, 0, 239847 },
                45737.7547,
                true
            },
            new object[] {
                "file 4",
                "98%HA",
                new List<double> { 0, 0 },
                45737.7547,
                false
            }
        };

        /// <summary>
        /// Test to check that stoichiometry is being calculated correctly in the for a peptide in the peptide:peptide case
        /// </summary>
        /// <param name="file">the filethe intensity data came from</param>
        /// <param name="group">the group the intensity is associated with</param>
        /// <param name="intensities">the intensities values for the numerator</param>
        /// <param name="intensity2">the intensity value of the denominator</param>
        /// <param name="useful">correct bool to compare to usefulStoichiometry</param>
        [Test]
        [TestCaseSource("_testPTMPeptidePeptideCase")]
        public void Stoichiometry_PTMPeptidePeptideCase_Pass(string file, string group, List<double> intensities,
            double intensity2, bool useful)
        {
            List<Intensity> Testintensities = new List<Intensity>();
            foreach (double d in intensities)
            {
                Testintensities.Add(new Intensity(file, group, d));
            }
            Intensity intensityTwo = new Intensity(file, group, intensity2);

            Stoichiometry TestStoichiometry = new Stoichiometry(Testintensities, intensityTwo);
            Assert.AreEqual(intensities.Sum() / intensity2, TestStoichiometry.StoichiometryVal);
            Assert.AreEqual(useful, TestStoichiometry.usefulStoichiometry);
        }

    }
}
