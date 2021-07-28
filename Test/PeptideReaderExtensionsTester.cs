using NUnit.Framework;
using System;
using System.Collections.Generic;
using PTMStoichiometry;
using System.Linq;
using System.IO;
using System.Text;

namespace Test
{
    [TestFixture]
    class PeptideReaderExtensionsTester
    {

        /// <summary>
        /// Test that checks that ReadTsv returns the correct number of Peptides
        /// </summary>
        /// <param name="filePath">file path to the file with the peptide data</param>
        /// <param name="groupPath">file path to the file with the group data</param>
        /// <param name="reqNumPepMeasurements">min num of peptide intensities that must be observed (non zero -> MS or MSMS detection)</param>
        /// <param name="peptideCount">the number of Peptides in the file</param>
        [Test]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-1000PeptidesProteinAlphabetized.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt", 3, 1000)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-5000PeptidesProteinAlphabetized.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt", 3, 5000)]
        [TestCase(@"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126-2021-07-07-08-59-09-AllQuantifiedPeptides-10000PeptidesProteinAlphabetized.txt",
            @"C:\Users\KAP\source\repos\PTMStoichiometry_master\Test\TestData\MSV000086126GlobalGroups.txt", 3, 10000)]
        public void PeptideReader_ReadTsv_Pass(string filePath, string groupPath, int reqNumPepMeasurements, int peptideCount)
        {
            //check throws error if file references no good
            //var noFilePeptideReaderTest = PeptideReader.ReadTsv("", filePathGroupsOnePeptideTest);
            //Assert.Throws(, onePeptideReaderTest.Count);

            List<Peptide> pepsInFile = PeptideReader.ReadTsv(filePath, groupPath, reqNumPepMeasurements, 5, "FlashLFQ");
            Assert.AreEqual(peptideCount, pepsInFile.Count());

            //test group dictonaries & lists!
        }

        

    }
}
