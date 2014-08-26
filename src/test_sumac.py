#!/usr/bin/env python


import unittest


class SumacTest(unittest.TestCase):
    """
    Tests for SUMAC
    """
    
    def test_version(self):
        import sys
        version = str(sys.version_info[0]) + "." + str(sys.version_info[1])
        self.assertEqual(version, "2.7", "SUMAC was written for Python 2.7. You are running Python " + version)

    
    
    def test_mafft_alignment(self):
        import os
        from alignments import Alignments
        alignments = Alignments(["test.fasta"], "unaligned")
        self.assertTrue(os.path.exists("alignments/test.fasta"))



    def setup_supermatrix(self):
        """
        Sets up supermatrix for some tests
        """
        from supermatrix import Supermatrix, Otu
        
        # make taxa
        alpha = Otu("alpha")
        beta = Otu("beta")
        pi = Otu("pi")
        omega = Otu("omega")
        gamma = Otu("gamma")

        # add sequence data to taxa, all that matters is the sequence length
        alpha_seq_lengths = [10, 0, 10, 0]
        beta_seq_lengths  = [5, 10, 0, 0]
        pi_seq_lengths    = [8, 6, 0, 4]
        omega_seq_lengths = [0, 0, 8, 9]
        gamma_seq_lengths = [0, 10, 0, 0]
        for length in alpha_seq_lengths:
            alpha.update("-", "x", length)
        for length in beta_seq_lengths:
            beta.update("-", "x", length)
        for length in pi_seq_lengths:
            pi.update("-", "x", length)
        for length in omega_seq_lengths:
            omega.update("-", "x", length)
        for length in gamma_seq_lengths:
            gamma.update("-", "x", length)
        
        sm = Supermatrix()
        sm.otus = {"alpha": alpha, "beta": beta, "pi": pi, "omega": omega, "gamma": gamma}
        sm.get_PD()
        return sm


        
    def test_pd_calculation(self):
        # 10 total triplets for 5 taxa, 2 of them are decisive, therefore pd=0.2
        sm = self.setup_supermatrix()
        self.assertEqual(sm.get_PD(), 0.2)



    def test_supermatrix_data_figure(self):
        import os
        sm = self.setup_supermatrix()
        sm.make_sequence_data_figure()
        self.assertTrue(os.path.exists("./sequence_data_plot.pdf"))



    def test_supermatrix_decisiveness_figure(self):
        import os
        sm = self.setup_supermatrix()
        sm.make_sequence_decisiveness_figure()
        self.assertTrue(os.path.exists("./sequence_decisiveness_plot_all_scores.pdf"))



    def test_supermatrix_decisiveness_csv(self):
        import os
        sm = self.setup_supermatrix()
        sm.make_decisiveness_csv()
        self.assertTrue(os.path.exists("./missing_sequence_decisiveness.csv"))


    #if verbose:
    #    print stuff

if __name__ == '__main__':
    unittest.main()
